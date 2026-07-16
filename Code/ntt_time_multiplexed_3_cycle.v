`timescale 1ns / 1ps
////////////////////////////////////////////////////////////////////////////////
// Company:     Indian Institute of Technology, Gandhinagar
// Engineer:    Pranay Arvind Patil 24110252
//
// Module Name: time_muxed_butterfly
// Description: 3-phase time-multiplexed NTT convolution over GF(257).
//              Computes c = INTT( NTT(a) o NTT(b) ) using a single shared
//              unified_transform_core over 3 clock cycles.
//
//              Phase 0 (cycle 1): NTT(a)   latch into ntt_a_reg
//              Phase 1 (cycle 2): NTT(b)   latch into ntt_b_reg
//              Phase 2 (cycle 3): INTT(dot_product)  latch into c
//
//              `done` pulses high for one cycle when c is valid.
////////////////////////////////////////////////////////////////////////////////
module butterfly #(
    parameter integer NUM_BUSES = 32,
    parameter integer BUS_WIDTH = 9
)(
    input  wire                               clk,
    input  wire                               rst,
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] a,
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] b,
    output reg  [(NUM_BUSES * BUS_WIDTH)-1:0] c,
    output reg                                done
);
    localparam integer LOG_NUM_BUSES = $clog2(NUM_BUSES);
    // Phase encoding
    localparam [1:0] PHASE_NTT_A = 2'd0;
    localparam [1:0] PHASE_NTT_B = 2'd1;
    localparam [1:0] PHASE_INTT  = 2'd2;
    reg [1:0] phase;
    // Storage for NTT outputs
    reg [(NUM_BUSES * BUS_WIDTH)-1:0] ntt_a_reg;
    reg [(NUM_BUSES * BUS_WIDTH)-1:0] ntt_b_reg;
    // Core control
    reg  [(NUM_BUSES * BUS_WIDTH)-1:0] core_in;
    reg                                core_is_inverse;
    wire [(NUM_BUSES * BUS_WIDTH)-1:0] core_out;
    // Dot product result (combinational, always driven from latched NTT outputs)
    wire [(NUM_BUSES * BUS_WIDTH)-1:0] dot_prod;
    // -------------------------------------------------------------------------
    // Shared transform core - one instance, used for all three phases
    // -------------------------------------------------------------------------
    unified_transform_core #(
        .NUM_BUSES (NUM_BUSES),
        .BUS_WIDTH (BUS_WIDTH)
    ) core_inst (
        .data_in    (core_in),
        .is_inverse (core_is_inverse),
        .data_out   (core_out)
    );
    // -------------------------------------------------------------------------
    // Dot product - combinationally driven from the two NTT result registers.
    // Includes N^{-1} * R correction so INTT output is in plain integer domain.
    // -------------------------------------------------------------------------
    dot_product #(
        .NUM_BUSES     (NUM_BUSES),
        .BUS_WIDTH     (BUS_WIDTH),
        .NUM_BUSES_LOG (LOG_NUM_BUSES)
    ) dot_product_inst (
        .a (ntt_a_reg),
        .b (ntt_b_reg),
        .c (dot_prod)
    );
    // -------------------------------------------------------------------------
    // Combinational input mux: feed the correct data into the core each phase
    // -------------------------------------------------------------------------
    always @(*) begin
        case (phase)
            PHASE_NTT_A: begin
                core_in         = a;
                core_is_inverse = 1'b0;
            end
            PHASE_NTT_B: begin
                core_in         = b;
                core_is_inverse = 1'b0;
            end
            PHASE_INTT: begin
                core_in         = dot_prod;  // INTT of the pointwise product
                core_is_inverse = 1'b1;
            end
            default: begin
                core_in         = a;
                core_is_inverse = 1'b0;
            end
        endcase
    end
    // -------------------------------------------------------------------------
    // FSM: latch core_out and advance phase on each clock edge
    // -------------------------------------------------------------------------
    always @(posedge clk or posedge rst) begin
        if (rst) begin
            phase     <= PHASE_NTT_A;
            done      <= 1'b0;
            ntt_a_reg <= {(NUM_BUSES * BUS_WIDTH){1'b0}};
            ntt_b_reg <= {(NUM_BUSES * BUS_WIDTH){1'b0}};
            c         <= {(NUM_BUSES * BUS_WIDTH){1'b0}};
        end else begin
            done <= 1'b0;   // default: deassert every cycle
            case (phase)
                PHASE_NTT_A: begin
                    ntt_a_reg <= core_out;      // latch NTT(a)
                    phase     <= PHASE_NTT_B;
                end
                PHASE_NTT_B: begin
                    ntt_b_reg <= core_out;      // latch NTT(b)
                    phase     <= PHASE_INTT;
                end
                PHASE_INTT: begin
                    c         <= core_out;      // latch final result
                    done      <= 1'b1;          // pulse done for one cycle
                    phase     <= PHASE_NTT_A;   // ready for next computation
                end
                default: phase <= PHASE_NTT_A;
            endcase
        end
    end
endmodule
////////////////////////////////////////////////////////////////////////////////
// Module Name: unified_transform_core
// Description: Shared combinational NTT/INTT core.
//              Uses Cooley-Tukey DIT for BOTH directions so the same physical
//              cooley_butterfly gates handle forward and inverse transforms.
//              `is_inverse` selects inverse twiddles at each butterfly.
//
// Data flow (fully combinational, zero latency):
//   data_in  bit_reverse_16  layer[1]  layer[2]  ...  layer[LOG]  data_out
////////////////////////////////////////////////////////////////////////////////
module unified_transform_core #(
    parameter integer NUM_BUSES = 32,
    parameter integer BUS_WIDTH = 16
)(
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] data_in,
    input  wire                               is_inverse,
    output wire [(NUM_BUSES * BUS_WIDTH)-1:0] data_out
);
    localparam integer LOG_NUM_BUSES = $clog2(NUM_BUSES);
    wire [(NUM_BUSES * BUS_WIDTH)-1:0] stage [0:LOG_NUM_BUSES];
    // Step 1: Bit-reversal - always applied for DIT regardless of direction
    bit_reverse_16 #(
        .NUM_BUSES (NUM_BUSES),
        .BUS_WIDTH (BUS_WIDTH)
    ) bit_rev_inst (
        .a (data_in),
        .c (stage[0])
    );
    // Step 2: Chain LOG_NUM_BUSES unified layers
    genvar L;
    generate
        for (L = 1; L <= LOG_NUM_BUSES; L = L + 1) begin : layers
            unified_cooley_layer #(
                .NUM_BUSES    (NUM_BUSES),
                .BUS_WIDTH    (BUS_WIDTH),
                .log_BUS_SIZE (LOG_NUM_BUSES),
                .layer        (L)
            ) layer_inst (
                .data_in    (stage[L-1]),
                .is_inverse (is_inverse),
                .data_out   (stage[L])
            );
        end
    endgenerate
    assign data_out = stage[LOG_NUM_BUSES];
endmodule
////////////////////////////////////////////////////////////////////////////////
// Module Name: unified_cooley_layer
// Description: One CT-DIT butterfly stage shared between forward NTT and
//              inverse INTT via a runtime twiddle MUX.
//
// SINGLE TABLE - HOW INVERSE TWIDDLES ARE DERIVED:
//   ROOTS[p] stores GEN^(256-p) * R mod Q  (p = 1..256, MSB = p=256)
//   where GEN=3, R=512, Q=257.
//
//   For butterfly exponent k (= TWIDDLE_EXP):
//     Forward twiddle  = GEN^k         P_FWD = 256 - k  (existing formula)
//     Inverse twiddle  = GEN^{-k}
//                      = GEN^{256-k}  (group order = 256)
//                       P_INV = 256 - (256-k) = k
//
//   Both indices are elaboration-time constants  zero runtime ROM access.
//   The `is_inverse` MUX costs BUS_WIDTH LUT2 cells per butterfly (9 LUTs
//   for BUS_WIDTH=9), negligible vs the butterfly multiply (~57 LUTs).
//
// LAYER PARAMETER MUST BE >= 1.
////////////////////////////////////////////////////////////////////////////////
module unified_cooley_layer #(
    parameter integer NUM_BUSES    = 8,
    parameter integer BUS_WIDTH    = 16,
    parameter integer log_BUS_SIZE = 3,
    parameter integer layer        = 1,
    parameter [4095:0] ROOTS = {
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
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] data_in,
    input  wire                               is_inverse,
    output wire [(NUM_BUSES * BUS_WIDTH)-1:0] data_out
);
    localparam integer STRIDE = 1 << layer;
    genvar start, j;
    generate
        for (start = 0; start < NUM_BUSES; start = start + STRIDE) begin : sub_ntt
            for (j = start; j < start + (STRIDE >> 1); j = j + 1) begin : butterfly_pos
                // -------------------------------------------------------------
                // Elaboration-time twiddle exponent for this butterfly position.
                // TWIDDLE_EXP = (j - start) * 2^(8 - layer)
                // When j == start: TWIDDLE_EXP = 0  twiddle = ω^0 = 1
                // -------------------------------------------------------------
                localparam integer TWIDDLE_EXP =
                    (j - start)
                    * (1 << (log_BUS_SIZE - layer))
                    * (1 << (8 - log_BUS_SIZE));
                // -------------------------------------------------------------
                // Forward index:  P_FWD = 256 - TWIDDLE_EXP
                //   Selects GEN^TWIDDLE_EXP * R from ROOTS table.
                // Inverse index:  P_INV = TWIDDLE_EXP  (or 256 when exp=0)
                //   Selects GEN^{-TWIDDLE_EXP} * R from same ROOTS table
                //   because GEN^{256-k} = GEN^{-k} in a group of order 256.
                // -------------------------------------------------------------
                localparam integer P_FWD = 256 - TWIDDLE_EXP;
                localparam integer P_INV = (TWIDDLE_EXP == 0) ? 256 : TWIDDLE_EXP;
                // Compile-time twiddle literals - no runtime ROM access
                localparam [BUS_WIDTH-1:0] w_fwd = ROOTS[16 * P_FWD - 1 -: 16];
                localparam [BUS_WIDTH-1:0] w_inv = ROOTS[16 * P_INV - 1 -: 16];
                // Runtime 2:1 MUX: BUS_WIDTH LUT2 cells, negligible overhead
                wire [BUS_WIDTH-1:0] w_sel = is_inverse ? w_inv : w_fwd;
                // -------------------------------------------------------------
                // Butterfly instantiation.
                // Bus convention: element i at [(NUM_BUSES-i)*BUS_WIDTH-1 -: BUS_WIDTH]
                // Upper pair: j,  Lower pair: j + STRIDE/2
                // -------------------------------------------------------------
                cooley_butterfly #(
                    .DATA_WIDTH (BUS_WIDTH)
                ) but_inst (
                    .a  (data_in [(NUM_BUSES - j)                   * BUS_WIDTH - 1 -: BUS_WIDTH]),
                    .b  (data_in [(NUM_BUSES - (j + (STRIDE >> 1))) * BUS_WIDTH - 1 -: BUS_WIDTH]),
                    .w  (w_sel),
                    .c0 (data_out[(NUM_BUSES - j)                   * BUS_WIDTH - 1 -: BUS_WIDTH]),
                    .c1 (data_out[(NUM_BUSES - (j + (STRIDE >> 1))) * BUS_WIDTH - 1 -: BUS_WIDTH])
                );
            end
        end
    endgenerate
endmodule
////////////////////////////////////////////////////////////////////////////////
// Module Name: cooley_butterfly     UNCHANGED
// Description: c0 = a + w*b mod Q,  c1 = a - w*b mod Q
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
    mont_mult #(
        .DATA_WIDTH (DATA_WIDTH)
    ) mult_inst (
        .a_mont (b),
        .b_mont (w),
        .c_mont (t)
    );
    wire [DATA_WIDTH:0] sum       = {1'b0, a} + {1'b0, t};
    wire [DATA_WIDTH:0] sum_sub_q = sum - {1'b0, Q};
    assign c0 = sum_sub_q[DATA_WIDTH] ? sum[DATA_WIDTH-1:0] : sum_sub_q[DATA_WIDTH-1:0];
    wire [DATA_WIDTH:0]   diff      = {1'b0, a} - {1'b0, t};
    wire [DATA_WIDTH-1:0] diff_add_q = diff[DATA_WIDTH-1:0] + Q;
    assign c1 = diff[DATA_WIDTH] ? diff_add_q : diff[DATA_WIDTH-1:0];
endmodule
////////////////////////////////////////////////////////////////////////////////
// Module Name: dot_product     UNCHANGED (CORR correction already included)
// Description: Pointwise Montgomery multiply + N^{-1} * R scaling.
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
    localparam integer CORR = (NINV * 512) % Q;
    genvar i;
    generate
        for (i = 0; i < NUM_BUSES; i = i + 1) begin : dot
            wire [BUS_WIDTH-1:0] raw;
            mont_mult #(
                .DATA_WIDTH (BUS_WIDTH)
            ) mult_inst (
                .a_mont (a[i*BUS_WIDTH +: BUS_WIDTH]),
                .b_mont (b[i*BUS_WIDTH +: BUS_WIDTH]),
                .c_mont (raw)
            );
            // Multiply by CORR (compile-time constant  shift-add tree)
            wire [12:0] prod = raw * CORR;
            // Barrett-style mod 257: (lo8) - (hi5), then add 257 if negative
            wire [9:0] diff     = {2'b0, prod[7:0]} - {5'b0, prod[12:8]};
            wire [9:0] diff_mod = diff[9] ? (diff + 10'd257) : diff;
            assign c[i*BUS_WIDTH +: BUS_WIDTH] = diff_mod;
        end
    endgenerate
endmodule
////////////////////////////////////////////////////////////////////////////////
//         DON'T TOUCH - Montgomery modular arithmetic primitives
////////////////////////////////////////////////////////////////////////////////
module mont_redc #(
    parameter integer DATA_WIDTH  = 16,
    parameter integer T_WIDTH     = 2 * DATA_WIDTH,
    parameter         Q           = 16'd257,
    parameter         Q_INV_R_NEG = 16'd255,
    parameter integer RSH         = 9
)(
    input  wire [T_WIDTH-1:0]    T,
    output wire [DATA_WIDTH-1:0] t_redc
);
    wire [RSH-1:0] m;
    generate
        if (Q == 16'd257 && Q_INV_R_NEG == 16'd255 && RSH == 9) begin : gen_opt_m
            assign m = {T[0], 8'h00} - T[RSH-1:0];
        end else begin : gen_generic_m
            assign m = T[RSH-1:0] * Q_INV_R_NEG[RSH-1:0];
        end
    endgenerate
    wire [DATA_WIDTH+RSH-1:0] P;
    generate
        if (Q == 16'd257 && RSH == 9) begin : gen_opt_P
            assign P = {1'b0, m, 8'h00} + m;
        end else begin : gen_generic_P
            assign P = m * Q;
        end
    endgenerate
    wire [DATA_WIDTH:0] t;
    generate
        if (Q == 16'd257 && RSH == 9) begin : gen_opt_t
            assign t = T[DATA_WIDTH+RSH-1:RSH] + P[DATA_WIDTH+RSH-1:RSH]
                       + (|T[RSH-1:0]);
        end else begin : gen_generic_t
            wire [T_WIDTH:0] TP = T + {{(T_WIDTH-DATA_WIDTH-RSH+1){1'b0}}, P};
            assign t = TP[DATA_WIDTH+RSH:RSH];
        end
    endgenerate
    wire [DATA_WIDTH:0] t_sub_q = t - {1'b0, Q[DATA_WIDTH-1:0]};
    assign t_redc = t_sub_q[DATA_WIDTH] ? t[DATA_WIDTH-1:0] : t_sub_q[DATA_WIDTH-1:0];
endmodule
module mont_mult #(
    parameter integer DATA_WIDTH  = 16,
    parameter integer T_WIDTH     = 2 * DATA_WIDTH,
    parameter         Q           = 16'd257,
    parameter         Q_INV_R_NEG = 16'd255,
    parameter integer RSH         = 9
)(
    input  wire [DATA_WIDTH-1:0] a_mont,
    input  wire [DATA_WIDTH-1:0] b_mont,
    output wire [DATA_WIDTH-1:0] c_mont
);
    (* use_dsp = "no" *) wire [T_WIDTH-1:0] T;
    generate
        if (Q == 16'd257 && RSH == 9) begin : gen_mult_opt
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
// Module Name: bit_reverse_16     UNCHANGED
////////////////////////////////////////////////////////////////////////////////
module bit_reverse_16 #(
    parameter NUM_BUSES = 8,
    parameter BUS_WIDTH = 16
)(
    input  wire [(NUM_BUSES * BUS_WIDTH)-1:0] a,
    output wire [(NUM_BUSES * BUS_WIDTH)-1:0] c
);
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