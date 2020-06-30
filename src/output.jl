using Printf

function output_result(stepnum,Qbase,cellxmax,cellymax,specific_heat_ratio,out_file_front,out_ext,out_dir)
    
    stepnum=string(stepnum)
    while length(stepnum) < 6
        stepnum="0"*stepnum
    end
    
    fff=out_dir*"/"*out_file_front*stepnum*out_ext
    open(fff,"w") do f
        write(f,"result:rho[kg/m^3], u[m/s], v[m/s], p[Pa], T[K]\n")
        for i in 2:cellxmax-1
            for j in 2:cellymax-1
                a1 = @sprintf("%8.8f", Qbase[i,j,1])
                a2 = @sprintf("%8.8f", Qbase[i,j,2])
                a3 = @sprintf("%8.8f", Qbase[i,j,3])
                a4 = @sprintf("%8.8f", Qbase[i,j,4])
                T = Qbase[i,j,4]*specific_heat_ratio/Qbase[i,j,1]
                a5 = @sprintf("%8.8f", T)
                write(f, a1*" "*a2*" "*a3*" "*a4*" "*a5*"\n")
            end
        end
    end
    println("\nwrite "*fff)
end
