`timescale 1ns / 1ps
////////////////////////////////////////////////////////////////////////////////
// Company:     Indian Institute of Technology, Gandhinagar
// Engineer:    Pranay Arvind Patil 24110252
//
// Module Name: butterfly
// Description: NTT-based negacyclic convolution over GF(257).
//              Computes c = INTT( NTT(a) o NTT(b) ) where o is pointwise
//              Montgomery multiplication, using Cooley-Tukey forward layers
//              and Gentleman-Sande inverse layers.
////////////////////////////////////////////////////////////////////////////////
module butterfly #(
    parameter integer NUM_BUSES = 32,   // Number of NTT lanes (must be a power of 2)
    parameter integer BUS_WIDTH = 9   // Bit-width of each lane
)(
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] a,
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] b,
    output wire [(NUM_BUSES * BUS_WIDTH)-1:0] c
);
    // log2(NUM_BUSES) - derived, not user-settable
    localparam integer LOG_NUM_BUSES = $clog2(NUM_BUSES);
    // -------------------------------------------------------------------------
    // Bit-reverse input permutation
    // -------------------------------------------------------------------------
    wire [(NUM_BUSES * BUS_WIDTH)-1:0] a_bit_reversed;
    wire [(NUM_BUSES * BUS_WIDTH)-1:0] b_bit_reversed;
    bit_reverse_16 #(
        .NUM_BUSES (NUM_BUSES),
        .BUS_WIDTH (BUS_WIDTH)
    ) reverser_inst1 (
        .a (a),
        .c (a_bit_reversed)
    );
    bit_reverse_16 #(
        .NUM_BUSES (NUM_BUSES),
        .BUS_WIDTH (BUS_WIDTH)
    ) reverser_inst2 (
        .a (b),
        .c (b_bit_reversed)
    );
    
    ///////// Forward cooley tukey ///////////////////
    wire [(NUM_BUSES * BUS_WIDTH)-1:0] stage_a [0:LOG_NUM_BUSES];
    wire [(NUM_BUSES * BUS_WIDTH)-1:0] stage_b [0:LOG_NUM_BUSES];
    assign stage_a[0] = a_bit_reversed;
    assign stage_b[0] = b_bit_reversed;
    genvar L;
    generate
        for (L = 1; L <= LOG_NUM_BUSES; L = L + 1) begin : layers
            dual_cooley_butterfly_layer #(
                .NUM_BUSES   (NUM_BUSES),
                .BUS_WIDTH   (BUS_WIDTH),
                .log_BUS_SIZE(LOG_NUM_BUSES),
                .layer       (L)
            ) layer_inst (
                .a_a (stage_a[L-1]),
                .a_b (stage_b[L-1]),
                .b_a (stage_a[L]),
                .b_b (stage_b[L])
            );
        end
    endgenerate
    // -------------------------------------------------------------------------
    // Pointwise dot product in Montgomery domain
    // -------------------------------------------------------------------------
    wire [(NUM_BUSES * BUS_WIDTH)-1:0] dot_prod;
    dot_product #(
        .NUM_BUSES (NUM_BUSES),
        .BUS_WIDTH (BUS_WIDTH),
        .NUM_BUSES_LOG (LOG_NUM_BUSES)
    ) dot_product_inst (
        .a (stage_a[LOG_NUM_BUSES]),
        .b (stage_b[LOG_NUM_BUSES]),
        .c (dot_prod)
    );
    localparam SCALE = ( 257 - (1 << (8 - LOG_NUM_BUSES))  ); // N^-1
    localparam CORR = (( 257 - (1 << (8 - LOG_NUM_BUSES))  ) * 512 ) %  257; // N^-1 * R     
    // -------------------------------------------------------------------------
    // Inverse NTT (Gentleman-Sande): stage LOG = dot product, stage 0 = output
    // -------------------------------------------------------------------------
    wire [(NUM_BUSES * BUS_WIDTH)-1:0] stage_gs [0:LOG_NUM_BUSES];
    assign stage_gs[LOG_NUM_BUSES] = dot_prod;
    genvar K;
    generate
        for (K = LOG_NUM_BUSES; K > 0; K = K - 1) begin : gs_layers
            gs_butterfly_layer #(
                .NUM_BUSES   (NUM_BUSES),
                .BUS_WIDTH   (BUS_WIDTH),
                .log_BUS_SIZE(LOG_NUM_BUSES),
                .layer       (K)
            ) layerGS_inst (
                .a (stage_gs[K]),
                .b (stage_gs[K-1])
            );
        end
    endgenerate
 
    // -------------------------------------------------------------------------
    // Output bit-reversal permutation
    // -------------------------------------------------------------------------
    bit_reverse_16 #(
        .NUM_BUSES (NUM_BUSES),
        .BUS_WIDTH (BUS_WIDTH)
    ) reverser_inst3 (
        .a (stage_gs[0]),
        .c (c)
    );
endmodule
////////////////////////////////////////////////////////////////////////////////
// Module Name: cooley_butterfly_layer
// Description: One stage of the Cooley-Tukey NTT.
//              LAYER PARAMETER MUST BE >= 1.
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Module Name: dual_cooley_butterfly_layer
// Description: One Cooley-Tukey NTT stage that processes TWO data streams (A
//              and B) simultaneously through a single generate loop.
//
//              KEY OPTIMIZATION vs two separate cooley_butterfly_layer instances:
//              - The twiddle factor `w` is a localparam computed ONCE per
//                butterfly position, then shared by both the stream-A and
//                stream-B cooley_butterfly units at that position.
//              - Both streams are synthesised in the same Vivado module context,
//                so constant propagation through mont_redc is identical for A
//                and B, eliminating the LUT asymmetry visible in the utilisation
//                report (e.g. layerA_inst=3320 vs layerB_inst=3606 at L=1).
//
// LAYER PARAMETER MUST BE >= 1.
////////////////////////////////////////////////////////////////////////////////
module dual_cooley_butterfly_layer #(
    parameter integer NUM_BUSES    = 8,
    parameter integer BUS_WIDTH    = 16,
    parameter integer log_BUS_SIZE = 3,
    parameter integer layer        = 2,
    parameter [4095:0] roots = {
        16'h00FF,16'h00FB,16'h00EF,16'h00CB,16'h005F,16'h001C,16'h0054,16'h00FC,
        16'h00F2,16'h00D4,16'h007A,16'h006D,16'h0046,16'h00D2,16'h0074,16'h005B,
        16'h0010,16'h0030,16'h0090,16'h00AF,16'h000B,16'h0021,16'h0063,16'h0028,
        16'h0078,16'h0067,16'h0034,16'h009C,16'h00D3,16'h0077,16'h0064,16'h002B,
        16'h0081,16'h0082,16'h0085,16'h008E,16'h00A9,16'h00FA,16'h00EC,16'h00C2,
        16'h0044,16'h00CC,16'h0062,16'h0025,16'h006F,16'h004C,16'h00E4,16'h00AA,
        16'h00FD,16'h00F5,16'h00DD,16'h0095,16'h00BE,16'h0038,16'h00A8,16'h00F7,
        16'h00E3,16'h00A7,16'h00F4,16'h00DA,16'h008C,16'h00A3,16'h00E8,16'h00B6,
        16'h0020,16'h0060,16'h001F,16'h005D,16'h0016,16'h0042,16'h00C6,16'h0050,
        16'h00F0,16'h00CE,16'h0068,16'h0037,16'h00A5,16'h00EE,16'h00C8,16'h0056,
        16'h0001,16'h0003,16'h0009,16'h001B,16'h0051,16'h00F3,16'h00D7,16'h0083,
        16'h0088,16'h0097,16'h00C4,16'h004A,16'h00DE,16'h0098,16'h00C7,16'h0053,
        16'h00F9,16'h00E9,16'h00B9,16'h0029,16'h007B,16'h0070,16'h004F,16'h00ED,
        16'h00C5,16'h004D,16'h00E7,16'h00B3,16'h0017,16'h0045,16'h00CF,16'h006B,
        16'h0040,16'h00C0,16'h003E,16'h00BA,16'h002C,16'h0084,16'h008B,16'h00A0,
        16'h00DF,16'h009B,16'h00D0,16'h006E,16'h0049,16'h00DB,16'h008F,16'h00AC,
        16'h0002,16'h0006,16'h0012,16'h0036,16'h00A2,16'h00E5,16'h00AD,16'h0005,
        16'h000F,16'h002D,16'h0087,16'h0094,16'h00BB,16'h002F,16'h008D,16'h00A6,
        16'h00F1,16'h00D1,16'h0071,16'h0052,16'h00F6,16'h00E0,16'h009E,16'h00D9,
        16'h0089,16'h009A,16'h00CD,16'h0065,16'h002E,16'h008A,16'h009D,16'h00D6,
        16'h0080,16'h007F,16'h007C,16'h0073,16'h0058,16'h0007,16'h0015,16'h003F,
        16'h00BD,16'h0035,16'h009F,16'h00DC,16'h0092,16'h00B5,16'h001D,16'h0057,
        16'h0004,16'h000C,16'h0024,16'h006C,16'h0043,16'h00C9,16'h0059,16'h000A,
        16'h001E,16'h005A,16'h000D,16'h0027,16'h0075,16'h005E,16'h0019,16'h004B,
        16'h00E1,16'h00A1,16'h00E2,16'h00A4,16'h00EB,16'h00BF,16'h003B,16'h00B1,
        16'h0011,16'h0033,16'h0099,16'h00CA,16'h005C,16'h0013,16'h0039,16'h00AB,
        16'h0100,16'h00FE,16'h00F8,16'h00E6,16'h00B0,16'h000E,16'h002A,16'h007E,
        16'h0079,16'h006A,16'h003D,16'h00B7,16'h0023,16'h0069,16'h003A,16'h00AE,
        16'h0008,16'h0018,16'h0048,16'h00D8,16'h0086,16'h0091,16'h00B2,16'h0014,
        16'h003C,16'h00B4,16'h001A,16'h004E,16'h00EA,16'h00BC,16'h0032,16'h0096,
        16'h00C1,16'h0041,16'h00C3,16'h0047,16'h00D5,16'h007D,16'h0076,16'h0061,
        16'h0022,16'h0066,16'h0031,16'h0093,16'h00B8,16'h0026,16'h0072,16'h0055
    }
)(
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] a_a,  // stream A input
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] a_b,  // stream B input
    output wire [(NUM_BUSES * BUS_WIDTH)-1:0] b_a,  // stream A output
    output wire [(NUM_BUSES * BUS_WIDTH)-1:0] b_b   // stream B output
);
    localparam integer stride = 1 << layer;
    genvar start, j;
    generate
        for (start = 0; start < NUM_BUSES; start = start + stride) begin : sub_ntt
            for (j = start; j < start + (stride >> 1); j = j + 1) begin : small_butterfly
                // Inside the generate loop, after computing w:
                localparam [BUS_WIDTH-1:0] w =
                    roots[16*(256 - ((j - start) * (1 << (log_BUS_SIZE - layer))) *
                          (1 << (8 - log_BUS_SIZE))) - 1 -: 16];
                
                // Add this guard for layer=1 where w is always 255:
                localparam [BUS_WIDTH-1:0] w_eff = (layer == 1) ? {BUS_WIDTH{1'b0}} | 8'd255 : w;
                // Stream A butterfly
                cooley_butterfly #(
                    .DATA_WIDTH (BUS_WIDTH)
                ) but_a (
                    .a  (a_a[(NUM_BUSES - j)                   * BUS_WIDTH - 1 -: BUS_WIDTH]),
                    .b  (a_a[(NUM_BUSES - (j + (stride >> 1))) * BUS_WIDTH - 1 -: BUS_WIDTH]),
                    .w  (w_eff),
                    .c0 (b_a[(NUM_BUSES - j)                   * BUS_WIDTH - 1 -: BUS_WIDTH]),
                    .c1 (b_a[(NUM_BUSES - (j + (stride >> 1))) * BUS_WIDTH - 1 -: BUS_WIDTH])
                );
                // Stream B butterfly - same w, different data
                cooley_butterfly #(
                    .DATA_WIDTH (BUS_WIDTH)
                ) but_b (
                    .a  (a_b[(NUM_BUSES - j)                   * BUS_WIDTH - 1 -: BUS_WIDTH]),
                    .b  (a_b[(NUM_BUSES - (j + (stride >> 1))) * BUS_WIDTH - 1 -: BUS_WIDTH]),
                    .w  (w_eff),
                    .c0 (b_b[(NUM_BUSES - j)                   * BUS_WIDTH - 1 -: BUS_WIDTH]),
                    .c1 (b_b[(NUM_BUSES - (j + (stride >> 1))) * BUS_WIDTH - 1 -: BUS_WIDTH])
                );
            end
        end
    endgenerate
endmodule
////////////////////////////////////////////////////////////////////////////////
// Module Name: gs_butterfly_layer
// Description: One stage of the Gentleman-Sande (inverse) NTT. Single-stream.
//              LAYER PARAMETER MUST BE >= 1.
////////////////////////////////////////////////////////////////////////////////
module gs_butterfly_layer #(
    parameter integer NUM_BUSES    = 8,
    parameter integer BUS_WIDTH    = 16,
    parameter integer log_BUS_SIZE = 3,
    parameter integer layer        = 2,
    parameter [4095:0] roots = {
        16'h00FF,16'h0055,16'h0072,16'h0026,16'h00B8,16'h0093,16'h0031,16'h0066,
        16'h0022,16'h0061,16'h0076,16'h007D,16'h00D5,16'h0047,16'h00C3,16'h0041,
        16'h00C1,16'h0096,16'h0032,16'h00BC,16'h00EA,16'h004E,16'h001A,16'h00B4,
        16'h003C,16'h0014,16'h00B2,16'h0091,16'h0086,16'h00D8,16'h0048,16'h0018,
        16'h0008,16'h00AE,16'h003A,16'h0069,16'h0023,16'h00B7,16'h003D,16'h006A,
        16'h0079,16'h007E,16'h002A,16'h000E,16'h00B0,16'h00E6,16'h00F8,16'h00FE,
        16'h0100,16'h00AB,16'h0039,16'h0013,16'h005C,16'h00CA,16'h0099,16'h0033,
        16'h0011,16'h00B1,16'h003B,16'h00BF,16'h00EB,16'h00A4,16'h00E2,16'h00A1,
        16'h00E1,16'h004B,16'h0019,16'h005E,16'h0075,16'h0027,16'h000D,16'h005A,
        16'h001E,16'h000A,16'h0059,16'h00C9,16'h0043,16'h006C,16'h0024,16'h000C,
        16'h0004,16'h0057,16'h001D,16'h00B5,16'h0092,16'h00DC,16'h009F,16'h0035,
        16'h00BD,16'h003F,16'h0015,16'h0007,16'h0058,16'h0073,16'h007C,16'h007F,
        16'h0080,16'h00D6,16'h009D,16'h008A,16'h002E,16'h0065,16'h00CD,16'h009A,
        16'h0089,16'h00D9,16'h009E,16'h00E0,16'h00F6,16'h0052,16'h0071,16'h00D1,
        16'h00F1,16'h00A6,16'h008D,16'h002F,16'h00BB,16'h0094,16'h0087,16'h002D,
        16'h000F,16'h0005,16'h00AD,16'h00E5,16'h00A2,16'h0036,16'h0012,16'h0006,
        16'h0002,16'h00AC,16'h008F,16'h00DB,16'h0049,16'h006E,16'h00D0,16'h009B,
        16'h00DF,16'h00A0,16'h008B,16'h0084,16'h002C,16'h00BA,16'h003E,16'h00C0,
        16'h0040,16'h006B,16'h00CF,16'h0045,16'h0017,16'h00B3,16'h00E7,16'h004D,
        16'h00C5,16'h00ED,16'h004F,16'h0070,16'h007B,16'h0029,16'h00B9,16'h00E9,
        16'h00F9,16'h0053,16'h00C7,16'h0098,16'h00DE,16'h004A,16'h00C4,16'h0097,
        16'h0088,16'h0083,16'h00D7,16'h00F3,16'h0051,16'h001B,16'h0009,16'h0003,
        16'h0001,16'h0056,16'h00C8,16'h00EE,16'h00A5,16'h0037,16'h0068,16'h00CE,
        16'h00F0,16'h0050,16'h00C6,16'h0042,16'h0016,16'h005D,16'h001F,16'h0060,
        16'h0020,16'h00B6,16'h00E8,16'h00A3,16'h008C,16'h00DA,16'h00F4,16'h00A7,
        16'h00E3,16'h00F7,16'h00A8,16'h0038,16'h00BE,16'h0095,16'h00DD,16'h00F5,
        16'h00FD,16'h00AA,16'h00E4,16'h004C,16'h006F,16'h0025,16'h0062,16'h00CC,
        16'h0044,16'h00C2,16'h00EC,16'h00FA,16'h00A9,16'h008E,16'h0085,16'h0082,
        16'h0081,16'h002B,16'h0064,16'h0077,16'h00D3,16'h009C,16'h0034,16'h0067,
        16'h0078,16'h0028,16'h0063,16'h0021,16'h000B,16'h00AF,16'h0090,16'h0030,
        16'h0010,16'h005B,16'h0074,16'h00D2,16'h0046,16'h006D,16'h007A,16'h00D4,
        16'h00F2,16'h00FC,16'h0054,16'h001C,16'h005F,16'h00CB,16'h00EF,16'h00FB
    }
)(
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] a,
    output wire [(NUM_BUSES * BUS_WIDTH)-1:0] b
);
    localparam integer stride = 1 << layer;
    genvar start, j;
    generate
        for (start = 0; start < NUM_BUSES; start = start + stride) begin : sub_ntt
            for (j = start; j < start + (stride >> 1); j = j + 1) begin : small_butterfly
                localparam integer twiddle_idx = 
                    256 - ((j - start) * (1 << (log_BUS_SIZE - layer))) *
                          (1 << (8 - log_BUS_SIZE));
                
                // roots[16*256 - 1 -: 16] = roots[4095:4080] = 8'hFF = 255
                // Vivado fails to fold this boundary slice as a constant.
                // Explicitly hardcode it when twiddle_idx hits 256.
                localparam [BUS_WIDTH-1:0] w_eff = 
                    (twiddle_idx == 256) ? {{(BUS_WIDTH-8){1'b0}}, 8'hFF} 
                                         : roots[16*twiddle_idx - 1 -: 16];
                gs_butterfly #(
                    .DATA_WIDTH (BUS_WIDTH)
                ) small_but_instance (
                    .a  (a[(NUM_BUSES - j)                   * BUS_WIDTH - 1 -: BUS_WIDTH]),
                    .b  (a[(NUM_BUSES - (j + (stride >> 1))) * BUS_WIDTH - 1 -: BUS_WIDTH]),
                    .w  (w_eff),
                    .c0 (b[(NUM_BUSES - j)                   * BUS_WIDTH - 1 -: BUS_WIDTH]),
                    .c1 (b[(NUM_BUSES - (j + (stride >> 1))) * BUS_WIDTH - 1 -: BUS_WIDTH])
                );
            end
        end
    endgenerate
endmodule
////////////////////////////////////////////////////////////////////////////////
// Module Name: cooley_butterfly
// Description: Cooley-Tukey butterfly: c0 = a + w*b mod Q, c1 = a - w*b mod Q
////////////////////////////////////////////////////////////////////////////////
module cooley_butterfly #(
    parameter integer DATA_WIDTH = 16
)(
    input  wire [DATA_WIDTH-1:0] a,
    input  wire [DATA_WIDTH-1:0] b,
    input  wire [DATA_WIDTH-1:0] w,
    output wire [DATA_WIDTH-1:0] c0,
    output wire [DATA_WIDTH-1:0] c1
);
    localparam [DATA_WIDTH-1:0] Q = 257;
    wire [DATA_WIDTH-1:0] t;
    // Montgomery Multiplication: t = w * b mod Q
    mont_mult #(
        .DATA_WIDTH (DATA_WIDTH)
    ) mult_inst (
        .a_mont (b),
        .b_mont (w),
        .c_mont (t)
    );
    // ==========================================
    // OPTIMIZATION 1: Modular Addition (c0)
    // ==========================================
    // Calculate sum and sum - Q concurrently.
    // We pad with 1'b0 to create a DATA_WIDTH+1 bit vector to catch the carry/borrow.
    wire [DATA_WIDTH:0] sum       = {1'b0, a} + {1'b0, t};
    wire [DATA_WIDTH:0] sum_sub_q = sum - {1'b0, Q};
    
    // The MSB of sum_sub_q acts as our comparator.
    // If sum_sub_q[DATA_WIDTH] is 1, a borrow occurred, meaning (a + t) < Q.
    // Therefore, we output the raw sum. Otherwise, it crossed Q, so we output sum_sub_q.
    assign c0 = sum_sub_q[DATA_WIDTH] ? sum[DATA_WIDTH-1:0] : sum_sub_q[DATA_WIDTH-1:0];
    // ==========================================
    // OPTIMIZATION 2: Modular Subtraction (c1)
    // ==========================================
    // Calculate a - t concurrently.
    wire [DATA_WIDTH:0] diff = {1'b0, a} - {1'b0, t};
    
    // If diff[DATA_WIDTH] is 1, a borrow occurred, meaning a < t. 
    // In modular arithmetic, negative results require adding Q.
    wire [DATA_WIDTH-1:0] diff_add_q = diff[DATA_WIDTH-1:0] + Q;
    
    // Mux selects between the standard difference or the Q-adjusted difference based on borrow.
    assign c1 = diff[DATA_WIDTH] ? diff_add_q : diff[DATA_WIDTH-1:0];
endmodule
///////////////////////////////////////////////////////////////////////////////
// Module Name: gs_butterfly
// Description: Gentleman-Sande butterfly: c0 = (a+b) mod Q, c1 = (a-b)*w mod Q
////////////////////////////////////////////////////////////////////////////////
module gs_butterfly #(
    parameter integer DATA_WIDTH = 16
)(
    input  wire [DATA_WIDTH-1:0] a,
    input  wire [DATA_WIDTH-1:0] b,
    input  wire [DATA_WIDTH-1:0] w,
    output wire [DATA_WIDTH-1:0] c0,
    output wire [DATA_WIDTH-1:0] c1
);
    localparam [DATA_WIDTH-1:0] Q = 257;
    // ==========================================
    // OPTIMIZATION 1: Modular Addition (c0)
    // ==========================================
    // Calculate sum and sum - Q concurrently.
    wire [DATA_WIDTH:0] sum       = {1'b0, a} + {1'b0, b};
    wire [DATA_WIDTH:0] sum_sub_q = sum - {1'b0, Q};
    
    // The MSB of sum_sub_q acts as our comparator.
    // If sum_sub_q[DATA_WIDTH] is 1, a borrow occurred, meaning (a + b) < Q.
    assign c0 = sum_sub_q[DATA_WIDTH] ? sum[DATA_WIDTH-1:0] : sum_sub_q[DATA_WIDTH-1:0];
    // ==========================================
    // OPTIMIZATION 2: Modular Subtraction (a - b)
    // ==========================================
    // Calculate a - b concurrently.
    wire [DATA_WIDTH:0] diff = {1'b0, a} - {1'b0, b};
    
    // If diff[DATA_WIDTH] is 1, a borrow occurred, meaning a < b. 
    wire [DATA_WIDTH-1:0] diff_add_q = diff[DATA_WIDTH-1:0] + Q;
    
    // Mux selects between the standard difference or the Q-adjusted difference.
    wire [DATA_WIDTH-1:0] diff_wire = diff[DATA_WIDTH] ? diff_add_q : diff[DATA_WIDTH-1:0];
    // ==========================================
    // Step 3: Montgomery Multiplication (c1)
    // ==========================================
    // c1 = (a - b) * w mod Q via Montgomery multiplication
    mont_mult #(
        .DATA_WIDTH (DATA_WIDTH)
    ) mult_inst (
        .a_mont (diff_wire),
        .b_mont (w),
        .c_mont (c1) // Optimization: Wire directly to the output port
    );
endmodule
////////////////////////////////////////////////////////////////////////////////
// Module Name: dot_product
// Description: Pointwise Montgomery multiplication of two NTT-domain vectors.
////////////////////////////////////////////////////////////////////////////////
module dot_product #(
    parameter integer NUM_BUSES     = 8,
    parameter integer BUS_WIDTH     = 16,
    parameter integer NUM_BUSES_LOG = 3
)(
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] a,
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] b,
    output wire [(NUM_BUSES * BUS_WIDTH)-1:0] c
);
    localparam integer Q    = 257;
    localparam integer NINV = Q - (1 << (8 - NUM_BUSES_LOG));
    localparam integer CORR = (NINV * 512) % Q;  // N^{-1} * R mod Q
    genvar i;
    generate
        for (i = 0; i < NUM_BUSES; i = i + 1) begin : dot
            wire [BUS_WIDTH-1:0] raw;
            // 1. Core Montgomery Multiplication
            mont_mult #(
                .DATA_WIDTH (BUS_WIDTH)
            ) mult_inst (
                .a_mont (a[i*BUS_WIDTH +: BUS_WIDTH]),
                .b_mont (b[i*BUS_WIDTH +: BUS_WIDTH]),
                .c_mont (raw)
            );
            // 2. Scale by CORR
            // Because CORR is a compile-time constant, Vivado natively infers
            // a highly optimized shift-add tree here (0 DSP slices).
            wire [12:0] prod = raw * CORR;
            // 3. O(1) Modulo 257 Reduction (The Fix)
            // X mod 257 = (X & 255) - (X >> 8)
            // We use 10 bits to safely handle the signed subtraction.
            wire [9:0] diff = {2'b0, prod[7:0]} - {5'b0, prod[12:8]};
            
            // If diff is negative (MSB is 1), add 257. Otherwise, output diff.
            wire [9:0] diff_mod = diff[9] ? (diff + 10'd257) : diff;
            // 4. Output Mapping
            // Implicit zero-extension to match BUS_WIDTH
            assign c[i*BUS_WIDTH +: BUS_WIDTH] = diff_mod;
        end
    endgenerate
endmodule
////////////////////////////////////////////////////////////////////////////////
//         DON'T TOUCH - Montgomery modular arithmetic primitives
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Module Name: mont_redc
// Description: Pure Montgomery REDC reduction.
//
// OPTIMISATIONS vs previous version:
//
// 1. NARROW ADDER  (biggest LUT saving in this module)
//    Previous: sum_full = T[T_WIDTH-1:RSH] + ...
//              T[T_WIDTH-1:RSH] = T[31:9] = 23-bit slice.
//              But for Q=257, T = 9×9 product  T[31:18] = 0 always.
//              Vivado propagated these zeros differently for stream A vs B,
//              causing the LUT asymmetry in the utilisation report.
//    Fixed:    t = T[DATA_WIDTH+RSH-1:RSH] + P[DATA_WIDTH+RSH-1:RSH] + carry
//              T[17:9] is the actual significant slice (9 bits, not 23).
//              This is a 9+9+1  10-bit adder vs the old 23+9+1  24-bit adder.
//              Both A and B now see exactly the same adder structure, fixing
//              the asymmetry.
//
// 2. TIGHTER CONDITIONAL SUBTRACTION
//    Previous used DATA_WIDTH+2 bits for t_sub. Now DATA_WIDTH+1 bits suffice
//    since t ≤ 2Q-1 < 2^(DATA_WIDTH).
////////////////////////////////////////////////////////////////////////////////
module mont_redc #(
    parameter integer DATA_WIDTH  = 16,
    parameter integer T_WIDTH = 2 * DATA_WIDTH,
    parameter         Q           = 16'd257,
    parameter         Q_INV_R_NEG = 16'd255,
    parameter integer RSH         = 9
)(
    input  wire [T_WIDTH-1:0]    T,
    output wire [DATA_WIDTH-1:0] t_redc
);
    // -------------------------------------------------------------------------
    // Step 1 - m = (T * Q_INV_R_NEG) mod R
    // Optimised for Q=257, Q_INV_R_NEG=255, RSH=9:
    //   255 mod 512 = 256 - 1, so m = T[8:0]*256 mod 512 - T[8:0]
    //                              = {T[0], 8'b0} - T[8:0]  (9-bit arithmetic)
    // -------------------------------------------------------------------------
    wire [RSH-1:0] m;
    generate
        if (Q == 16'd257 && Q_INV_R_NEG == 16'd255 && RSH == 9) begin : gen_opt_m
            assign m = {T[0], 8'h00} - T[RSH-1:0];
        end else begin : gen_generic_m
            assign m = T[RSH-1:0] * Q_INV_R_NEG[RSH-1:0];
        end
    endgenerate
    // -------------------------------------------------------------------------
    // Step 2 - P = m * Q
    // Optimised for Q=257, RSH=9:
    //   m * 257 = m * 256 + m = {m, 8'b0} + m  (shift-add, no multiplier)
    // -------------------------------------------------------------------------
    wire [DATA_WIDTH+RSH-1:0] P;
    generate
        if (Q == 16'd257 && RSH == 9) begin : gen_opt_P
            assign P = {1'b0, m, 8'h00} + m;
        end else begin : gen_generic_P
            assign P = m * Q;
        end
    endgenerate
    // -------------------------------------------------------------------------
    // Step 3 - t = (T + P) >> RSH
    //
    // CARRY CORRECTNESS for Q=257:
    //   By Montgomery's property, (T + P) mod R = 0 exactly, so the lower RSH
    //   bits of (T + P) are always zero. The carry out of those bits into
    //   bit RSH equals (T[RSH-1:0] != 0), i.e. (|T[RSH-1:0]).
    //   Proof: P[RSH-1:0] = (R - T[RSH-1:0]) mod R, so T[RSH-1:0] + P[RSH-1:0]
    //   is either 0 (when T[RSH-1:0]=0, carry=0) or R (carry=1).
    //
    // ADDER WIDTH FIX:
    //   Use T[DATA_WIDTH+RSH-1:RSH] (9 bits when DATA_WIDTH=9) instead of
    //   T[T_WIDTH-1:RSH] (23 bits). T[T_WIDTH-1:DATA_WIDTH+RSH] = 0 always
    //   for a 9×9 product, but Vivado cannot prove that from the wire widths
    //   alone, so it synthesises the full 23-bit adder and optimises each
    //   cooley_butterfly_layer instance differently  LUT asymmetry.
    //   Constraining to 9 bits removes this ambiguity entirely.
    // -------------------------------------------------------------------------
    wire [DATA_WIDTH:0] t;  // DATA_WIDTH+1 bits; t < 2Q so fits safely
    generate
        if (Q == 16'd257 && RSH == 9) begin : gen_opt_t
            assign t = T[DATA_WIDTH+RSH-1:RSH] + P[DATA_WIDTH+RSH-1:RSH]
                       + (|T[RSH-1:0]);
        end else begin : gen_generic_t
            // Generic: full-width addition then shift
            wire [T_WIDTH:0] TP = T + {{(T_WIDTH-DATA_WIDTH-RSH+1){1'b0}}, P};
            assign t = TP[DATA_WIDTH+RSH:RSH];
        end
    endgenerate
    // -------------------------------------------------------------------------
    // Step 4 - conditional subtraction: output t if t < Q, else t - Q
    // DATA_WIDTH+1 bits suffice (t < 2Q, so t - Q fits in DATA_WIDTH bits).
    // -------------------------------------------------------------------------
    wire [DATA_WIDTH:0] t_sub_q = t - {1'b0, Q[DATA_WIDTH-1:0]};
    assign t_redc = t_sub_q[DATA_WIDTH] ? t[DATA_WIDTH-1:0] : t_sub_q[DATA_WIDTH-1:0];
endmodule
////////////////////////////////////////////////////////////////////////////////
// Module Name: mont_mult
// Description: Montgomery multiplication.
//
// OPTIMISATION: (* use_dsp = "yes" *) attribute on the product wire.
// The ZCU104 (UltraScale+) has 2520 DSP58E2 blocks. Each DSP handles an
// 18×27 multiply natively. For Q=257 (9-bit operands) one DSP per mont_mult
// is sufficient, freeing the CARRY8 and LUT resources previously used.
// At ~370 mont_mult instances in a 32-point NTT this uses <15% of DSPs.
// Remove the attribute if you prefer LUT-only (e.g. for DSP-constrained builds).
////////////////////////////////////////////////////////////////////////////////
module mont_mult #(
    parameter integer DATA_WIDTH  = 16,
    parameter integer T_WIDTH = 2 * DATA_WIDTH,
    parameter         Q           = 16'd257,
    parameter         Q_INV_R_NEG = 16'd255,
    parameter integer RSH         = 9
)(
    input  wire [DATA_WIDTH-1:0] a_mont,
    input  wire [DATA_WIDTH-1:0] b_mont,
    output wire [DATA_WIDTH-1:0] c_mont
);
    // Hint to Vivado: map the multiply onto a DSP block.
    // For Q=257, RSH=9 the operands are 9 bits - one DSP58E2 per instance.
    (* use_dsp = "no" *) wire [T_WIDTH-1:0] T;
    //wire [T_WIDTH-1:0] T;
    generate
        if (Q == 16'd257 && RSH == 9) begin : gen_mult_opt
            // Explicitly 9×9: Vivado sees a narrow multiply and assigns one DSP.
            assign T = a_mont[RSH-1:0] * b_mont[RSH-1:0];   
        end else begin : gen_mult_generic
            assign T = a_mont * b_mont;
        end
    endgenerate
    mont_redc #(
        .DATA_WIDTH  (DATA_WIDTH),
        .T_WIDTH     (T_WIDTH),
        .Q           (Q),
        .Q_INV_R_NEG (Q_INV_R_NEG),
        .RSH         (RSH)
    ) redc_inst (
        .T      (T),
        .t_redc (c_mont)
    );
endmodule
////////////////////////////////////////////////////////////////////////////////
// Module Name: mont_decode
// Description: Montgomery domain  integer domain.
////////////////////////////////////////////////////////////////////////////////
module mont_decode #(
    parameter integer DATA_WIDTH  = 16,
    parameter integer T_WIDTH = 2 * DATA_WIDTH,
    parameter         Q           = 257,
    parameter         Q_INV_R_NEG = 255,
    parameter integer RSH         = 9
)(
    input  wire [DATA_WIDTH-1:0] x_mont,
    output wire [DATA_WIDTH-1:0] x_int
);
    wire [T_WIDTH-1:0] t_in = { {(T_WIDTH - DATA_WIDTH){1'b0}}, x_mont };
    mont_redc #(
        .DATA_WIDTH  (DATA_WIDTH),
        .T_WIDTH     (T_WIDTH),
        .Q           (Q),
        .Q_INV_R_NEG (Q_INV_R_NEG),
        .RSH         (RSH)
    ) redc_inst (
        .T     (t_in),
        .t_redc(x_int)
    );
endmodule
////////////////////////////////////////////////////////////////////////////////
// Module Name: mont_encode
// Description: Integer domain  Montgomery domain.
////////////////////////////////////////////////////////////////////////////////
module mont_encode #(
    parameter integer DATA_WIDTH  = 16,
    parameter integer T_WIDTH     = 2 * DATA_WIDTH,
    parameter         Q           = 257,
    parameter         Q_INV_R_NEG = 255,
    parameter integer RSH         = 9,
    parameter         R2_MOD_Q    = 4
)(
    input  wire [DATA_WIDTH-1:0] x_int,
    output wire [DATA_WIDTH-1:0] x_mont
);
    // R2_MOD_Q = 4 = 2^2, so x_int * 4 is a left shift by 2.
    // No multiplier needed - just wire reassignment.
    wire [T_WIDTH-1:0] T = {{(T_WIDTH-DATA_WIDTH-2){1'b0}}, x_int, 2'b00};
    mont_redc #(
        .DATA_WIDTH  (DATA_WIDTH),
        .T_WIDTH     (T_WIDTH),
        .Q           (Q),
        .Q_INV_R_NEG (Q_INV_R_NEG),
        .RSH         (RSH)
    ) redc_inst (
        .T      (T),
        .t_redc (x_mont)
    );
endmodule
////////////////////////////////////////////////////////////////////////////////
// Company:     Indian Institute of Technology, Gandhinagar
// Engineer:    Pranay Arvind Patil 24110252
//
// Module Name: bit_reverse_16
// Description: Bit-reversal permutation for NUM_BUSES lanes of BUS_WIDTH bits.
//              Lane index i is routed to its bit-reversed index rev(i).
////////////////////////////////////////////////////////////////////////////////
module bit_reverse_16 #(
    parameter NUM_BUSES = 8,
    parameter BUS_WIDTH = 16
)(
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] a,
    output wire [(NUM_BUSES * BUS_WIDTH)-1:0] c
);
    // Verilog-2001 equivalent of $clog2: returns ceil(log2(value))
    function integer clog2;
        input integer value;
        integer i;
        begin
            clog2 = 0;
            for (i = value - 1; i > 0; i = i >> 1)
                clog2 = clog2 + 1;
        end
    endfunction
    localparam ADDR_BITS = clog2(NUM_BUSES);
    // Returns the ADDR_BITS-wide bit-reversal of idx
    function integer bit_reverse;
        input integer idx;
        integer r, j;
        begin
            r = 0;
            for (j = 0; j < ADDR_BITS; j = j + 1)
                r = (r << 1) | ((idx >> j) & 1);
            bit_reverse = r;
        end
    endfunction
    genvar gv;
    generate
        for (gv = 0; gv < NUM_BUSES; gv = gv + 1) begin : bitrev
            localparam rev = bit_reverse(gv);
            assign c[gv * BUS_WIDTH +: BUS_WIDTH] = a[rev * BUS_WIDTH +: BUS_WIDTH];
        end
    endgenerate
endmodule