`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: IIT Gandhinagar 
// Engineer: Pranay Patil
// 
// Create Date: 03/04/2026 11:24:40 AM
// Design Name: ntt_sdf
// Module Name: trial_modules
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////



module trial_modules #(
    parameter integer NUM_BUSES = 32,   // Number of NTT lanes (must be a power of 2)
    parameter integer BUS_WIDTH = 9   // Bit-width of each lane
)(
    input  wire [BUS_WIDTH - 1: 0] a, // Streaming input: 1 word per clock, already bit-reversed
    input  wire [BUS_WIDTH - 1: 0] b, // Streaming input: 1 word per clock, already bit-reversed
    input  wire clk,
    input  wire reset,
    output wire [BUS_WIDTH - 1: 0] c  // Streaming output: 1 word per clock, normal order
);

    localparam integer LOG_NUM_BUSES = $clog2(NUM_BUSES);

    // ------------------------------------------------
    // Global Counter (Controls Streaming and Routing)
    // ------------------------------------------------
    reg [LOG_NUM_BUSES - 1: 0] global_counter;
    always @(posedge clk) begin
        if (reset) begin
            global_counter <= 0;
        end else begin
            global_counter <= global_counter + 1; 
        end
    end

    // ------------------------------------------------
    // The Pipelines
    // ------------------------------------------------
    wire [BUS_WIDTH - 1: 0] pipe_ct_a [0:LOG_NUM_BUSES];
    wire [BUS_WIDTH - 1: 0] pipe_ct_b [0:LOG_NUM_BUSES];
    
    // The streaming inputs plug directly into the very first stage of the pipeline!
    assign pipe_ct_a[0] = a;
    assign pipe_ct_b[0] = b;

    // ----------------------------------------------
    // Definition of the SDF unit CT (Dual Pipeline)
    // ----------------------------------------------
    genvar index0;
    generate
        for(index0 = 1; index0 <= LOG_NUM_BUSES; index0 = index0 + 1) begin: BU_unit_ct

            wire sel = global_counter[index0 - 1];

            reg [LOG_NUM_BUSES - 1: 0] bu_counter;
            always @(posedge clk) begin
                if (reset) bu_counter <= 0;
                else if (sel) bu_counter <= bu_counter + 1;
            end

            // ==========================================
            // DATAPATH A
            // ==========================================
            wire [BUS_WIDTH - 1: 0] stage_in_a  = pipe_ct_a[index0 - 1];
            wire [BUS_WIDTH - 1: 0] stage_out_a;            
            wire [BUS_WIDTH - 1: 0] reg_in_a;
            wire [BUS_WIDTH - 1: 0] reg_out_a;
            wire [BUS_WIDTH - 1: 0] BU_out0_a;
            wire [BUS_WIDTH - 1: 0] BU_out1_a;
            
            BU_CT #(.BUS_WIDTH(BUS_WIDTH), .LOG_NUM_BUSES(LOG_NUM_BUSES), .layer(index0)) BU_CT_inst_a (
                .a(reg_out_a),
                .b(stage_in_a),
                .counter(bu_counter), 
                .c0(BU_out0_a),
                .c1(BU_out1_a)
            );

            shift_register #(.WIDTH(BUS_WIDTH), .DEPTH(1 << (index0-1))) shift_reg_inst_a (
                .clk(clk),
                .in_data(reg_in_a),
                .out_data(reg_out_a)
            );

             assign stage_out_a = sel ? BU_out0_a : reg_out_a;
             assign reg_in_a    = sel ? BU_out1_a : stage_in_a;
             assign pipe_ct_a[index0] = stage_out_a;

            // ==========================================
            // DATAPATH B
            // ==========================================
            wire [BUS_WIDTH - 1: 0] stage_in_b  = pipe_ct_b[index0 - 1];
            wire [BUS_WIDTH - 1: 0] stage_out_b;            
            wire [BUS_WIDTH - 1: 0] reg_in_b;
            wire [BUS_WIDTH - 1: 0] reg_out_b;
            wire [BUS_WIDTH - 1: 0] BU_out0_b;
            wire [BUS_WIDTH - 1: 0] BU_out1_b;
            
            BU_CT #(.BUS_WIDTH(BUS_WIDTH), .LOG_NUM_BUSES(LOG_NUM_BUSES), .layer(index0)) BU_CT_inst_b (
                .a(reg_out_b),
                .b(stage_in_b),
                .counter(bu_counter), 
                .c0(BU_out0_b),
                .c1(BU_out1_b)
            );

            shift_register #(.WIDTH(BUS_WIDTH), .DEPTH(1 << (index0-1))) shift_reg_inst_b (
                .clk(clk),
                .in_data(reg_in_b),
                .out_data(reg_out_b)
            );

             assign stage_out_b = sel ? BU_out0_b : reg_out_b;
             assign reg_in_b    = sel ? BU_out1_b : stage_in_b;
             assign pipe_ct_b[index0] = stage_out_b;

        end
    endgenerate

    wire [BUS_WIDTH - 1: 0] final_ct_out_a = pipe_ct_a[LOG_NUM_BUSES];
    wire [BUS_WIDTH - 1: 0] final_ct_out_b = pipe_ct_b[LOG_NUM_BUSES];

    // -----------------------------------------------
    // Definition of the Pointwise Multiplier
    // -----------------------------------------------
    wire [BUS_WIDTH - 1: 0] dot_prod_out;

    dot_product #(
        .BUS_WIDTH(BUS_WIDTH),
        .NUM_BUSES_LOG(LOG_NUM_BUSES)
    ) pointwise_mult_inst (
        .a(final_ct_out_a),
        .b(final_ct_out_b),
        .c(dot_prod_out)
    );

    // -----------------------------------------------
    // Definition of SDF GS unit (Decimation-in-Frequency)
    // -----------------------------------------------
    wire [BUS_WIDTH - 1: 0] pipe_gs [0:LOG_NUM_BUSES];
    assign pipe_gs[0] = dot_prod_out; 

    genvar index1;
    generate
        for(index1 = 1; index1 <= LOG_NUM_BUSES; index1 = index1 + 1) begin: BU_unit_gs

            wire sel = global_counter[LOG_NUM_BUSES - index1];

            reg [LOG_NUM_BUSES - 1: 0] bu_counter_gs;
            always @(posedge clk) begin
                if (reset) bu_counter_gs <= 0;
                else if (sel) bu_counter_gs <= bu_counter_gs + 1;
            end

            wire [BUS_WIDTH - 1: 0] stage_in  = pipe_gs[index1 - 1];
            wire [BUS_WIDTH - 1: 0] stage_out;            
            wire [BUS_WIDTH - 1: 0] reg_in;
            wire [BUS_WIDTH - 1: 0] reg_out;
            wire [BUS_WIDTH - 1: 0] BU_out0;
            wire [BUS_WIDTH - 1: 0] BU_out1;
            
            BU_GS #(.BUS_WIDTH(BUS_WIDTH), .LOG_NUM_BUSES(LOG_NUM_BUSES), .layer(index1)) BU_GS_inst (
                .a(reg_out),
                .b(stage_in),
                .counter(bu_counter_gs), 
                .c0(BU_out0),
                .c1(BU_out1)
            );

            shift_register #(.WIDTH(BUS_WIDTH), .DEPTH(1 << (LOG_NUM_BUSES - index1))) shift_reg_inst (
                .clk(clk),
                .in_data(reg_in), 
                .out_data(reg_out)
            );

             assign stage_out = sel ? BU_out0 : reg_out;
             assign reg_in    = sel ? BU_out1 : stage_in;
             assign pipe_gs[index1] = stage_out;

        end
    endgenerate

    // The final result streams directly out of the last stage and into the output pin!
    assign c = pipe_gs[LOG_NUM_BUSES];

endmodule




module shift_register #(
    parameter integer WIDTH = 9,
    parameter integer DEPTH = 128 // Must be a power of 2
)(
    input  wire             clk,
    input  wire [WIDTH-1:0] in_data,
    output wire [WIDTH-1:0] out_data
);

    // Automatically calculate how many bits we need for the address pointer
    localparam PTR_WIDTH = $clog2(DEPTH);
    
    // Force Vivado to use BRAM for large delay lines
    (* ram_style = "block" *) reg [WIDTH-1:0] ram [0:DEPTH-1];
    
    reg [PTR_WIDTH-1:0] ptr = 0;
    reg [WIDTH-1:0]     out_reg = 0;

    always @(posedge clk) begin
        // Read the oldest data at the current pointer
        out_reg <= ram[ptr]; 
        
        // Overwrite it with the newest data
        ram[ptr] <= in_data; 
        
        // Move the pointer forward (it automatically wraps to 0 because of bit-width)
        ptr <= ptr + 1;
    end
    
    assign out_data = out_reg;

endmodule



module BU_CT #(
    parameter integer NUM_BUSES = 32,   // Number of NTT lanes (must be a power of 2)
    parameter integer BUS_WIDTH = 9,   // Bit-width of each lane
    parameter integer LOG_NUM_BUSES = 8,
    parameter integer CT_GS = 0,
    parameter [4:0] layer = 2,
    parameter Q = 257,
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
        16'h0010,16'h005B,16'h0074,16'h00D2,16'h0046,16'h006D,16'h007A,16'h00D4
     }

)(
    input wire [BUS_WIDTH - 1: 0] a,
    input wire [BUS_WIDTH - 1: 0] b,
    input wire [LOG_NUM_BUSES - 1: 0] counter,
    output wire [BUS_WIDTH - 1: 0] c0,
    output wire [BUS_WIDTH - 1: 0] c1
    
);
    localparam stride = 1 << (layer);

    wire [BUS_WIDTH-1:0] w_fw = roots[16*(256 - (    ( counter & ( (stride >> 1) - 1) ) * ( 1 << (LOG_NUM_BUSES - layer) ) * ( 1 <<  (8 - LOG_NUM_BUSES)  )   )  ) - 1 -: 16];

    // counter will dictate what root will be chosen, so layer and the counter are the inputs, in that in a single stride, we care about counter % stride is of interest to us because the roots repeat after
    // each stride. w^{i}_(n) = w^{1}_{n/i}; w^{i}_{n/(2^(3-layer))} is the root that we choose, where i belongs to [0, n/2); This is implemented in the indexing
    // also (1 << (8 - LOG_NUM_BUSES)) is multiplied to reduce the 256 roots to n roots. Like we skip the unncessary roots and only the important ones remain
    
    wire [BUS_WIDTH - 1: 0] t;
    mont_mult #(.DATA_WIDTH(BUS_WIDTH)) mult_inst (.a_mont( w_fw) , .b_mont(b), .c_mont(t));

    // we need to add another bit to the t and a wires to ensure that overflow doesn't happen and cause the comparision to evaluate to the wrong value.
    wire [BUS_WIDTH: 0] sum = {1'b0, a} + {1'b0, t};
    wire [BUS_WIDTH: 0] sum_sub_Q = sum - Q;
    assign c0 = sum_sub_Q[BUS_WIDTH] ? sum[BUS_WIDTH-1:0] : sum_sub_Q[BUS_WIDTH-1:0];


    wire [BUS_WIDTH: 0] diff = {1'b0, a} - {1'b0, t};
    wire [BUS_WIDTH - 1: 0] diff_add_q = diff[BUS_WIDTH - 1: 0] + Q;
    assign c1 = diff[BUS_WIDTH] ? diff_add_q : diff[BUS_WIDTH - 1: 0];


endmodule


module BU_GS #(
    parameter integer NUM_BUSES = 32,   // Number of NTT lanes (must be a power of 2)
    parameter integer BUS_WIDTH = 9,   // Bit-width of each lane
    parameter integer LOG_NUM_BUSES = 8,
    parameter integer CT_GS = 0,
    parameter [4:0] layer = 2,
    parameter Q = 257,
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
        16'h0010,16'h005B,16'h0074,16'h00D2,16'h0046,16'h006D,16'h007A,16'h00D4
     }

)(
    input wire [BUS_WIDTH - 1: 0] a,
    input wire [BUS_WIDTH - 1: 0] b,
    input wire [LOG_NUM_BUSES - 1: 0] counter,
    output wire [BUS_WIDTH - 1: 0] c0,
    output wire [BUS_WIDTH - 1: 0] c1
    
);  
    localparam stride = 1 << (LOG_NUM_BUSES - layer + 1);
    
    wire [7:0] index = ( counter & ( (stride >> 1) - 1) ) * ( 1 << (LOG_NUM_BUSES - layer) ) * ( 1 <<  (8 - LOG_NUM_BUSES)  ); 

    // FIXED: Added (index + 1) to prevent the -1 out-of-bounds crash when index is 0
    wire [BUS_WIDTH-1:0] w_in = roots[16 * (index + 1) - 1 -: 16];

    // counter will dictate what root will be chosen, so layer and the counter are the inputs, in that in a single stride, we care about counter % stride is of interest to us because the roots repeat after
    // each stride. w^{i}_(n) = w^{1}_{n/i}; w^{i}_{n/(2^(3 - layer))} is the root that we choose, where i belongs to [0, n/2); This is implemented in the indexing
    // also (1 << (8 - LOG_NUM_BUSES)) is multiplied to reduce the 256 roots to n roots. Like we skip the unncessary roots and only the important ones remain
    // and w^-1 = w^(N - 1)


    // we need to add another bit to the t and a wires to ensure that overflow doesn't happen and cause the comparision to evaluate to the wrong value.
    wire [BUS_WIDTH: 0] sum = {1'b0, a} + {1'b0, b};
    wire [BUS_WIDTH: 0] sum_sub_Q = sum - Q;
    assign c0 = sum_sub_Q[BUS_WIDTH] ? sum[BUS_WIDTH-1:0] : sum_sub_Q[BUS_WIDTH-1:0];


    wire [BUS_WIDTH: 0] diff = {1'b0, a} - {1'b0, b};
    wire [BUS_WIDTH - 1: 0] diff_add_q = diff[BUS_WIDTH - 1: 0] + Q;
    wire [BUS_WIDTH - 1: 0] t;
    assign t = diff[BUS_WIDTH] ? diff_add_q : diff[BUS_WIDTH - 1: 0];
    mont_mult #(.DATA_WIDTH( BUS_WIDTH )) mult_inst (.a_mont( w_in ) , .b_mont( t ), .c_mont( c1 ));



endmodule



////////////////////////////////////////////////////////////////////////////////
// Module Name: scalar_dot_product
// Description: Streaming Pointwise Montgomery mult with N^-1 scaling
////////////////////////////////////////////////////////////////////////////////
module dot_product #(
    parameter integer BUS_WIDTH     = 9,
    parameter integer NUM_BUSES_LOG = 5
)(
    input  wire [BUS_WIDTH-1:0] a,
    input  wire [BUS_WIDTH-1:0] b,
    output wire [BUS_WIDTH-1:0] c
);
    localparam integer Q    = 257;
    // Calculate N^-1 scaling factor based on the pipeline depth
    localparam integer NINV = Q - (1 << (8 - NUM_BUSES_LOG));
    localparam integer CORR = (NINV * 512) % Q;  

    wire [BUS_WIDTH-1:0] raw;

    // 1. Core Montgomery Multiplication
    mont_mult #(
        .DATA_WIDTH (BUS_WIDTH)
    ) mult_inst (
        .a_mont (a),
        .b_mont (b),
        .c_mont (raw)
    );

    // 2. Scale by CORR (0 DSP slices, just a shift-add tree)
    wire [12:0] prod = raw * CORR;

    // 3. O(1) Modulo 257 Reduction
    wire [9:0] diff = {2'b0, prod[7:0]} - {5'b0, prod[12:8]};
    wire [9:0] diff_mod = diff[9] ? (diff + 10'd257) : diff;

    // 4. Output Mapping
    assign c = diff_mod[BUS_WIDTH-1:0];

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