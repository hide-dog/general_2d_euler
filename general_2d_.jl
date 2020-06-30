# src path
# src_path="C:\\Users\\秀人\\Downloads\\一般二次元オイラー方程式\\General_euler\\src\\"
src_path="C:\\Users\\Koike\\Google ドライブ\\Programing言語\\programing_hobby\\一般二次元圧縮性オイラー方程式\\General_euler\\src\\"

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
src_main="main.jl"
include(src_path*src_main)


