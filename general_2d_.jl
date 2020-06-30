# src path
src_path="C:\\Users\\Koike\\Desktop\\git\\general_2d_euler\\src\\"

# main (変更しないこと)
src_read="read_grid.jl"
include(src_path*src_read)
src_read="read_para.jl"
include(src_path*src_read)
src_read="cal_time_step.jl"
include(src_path*src_read)
src_read="output.jl"
include(src_path*src_read)
src_read="setup.jl"
include(src_path*src_read)
src_read="boundary.jl"
include(src_path*src_read)
src_read="misc.jl"
include(src_path*src_read)
src_main="main.jl"
include(src_path*src_main)


