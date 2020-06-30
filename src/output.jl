function output_result(stepnum,Qbase,specific_heat_ratio,out_file_front,out_ext,out_dir) 
    stepnum=string(stepnum)

    while length(stepnum) < 6
        stepnum="0"*stepnum
    end
    
    fff=out_dir*"/"*out_file_front*stepnum*out_ext
    open(fff,"w") do f
        write(f,"result:rho[kg/m^3], u[m/s], v[m/s], p[Pa], T[K]\n")
        for i in 1:size(Qbase)[1]
            T=Qbase[i,4]*specific_heat_ratio/Qbase[i,1]
            write(f,string(i)*" "*string(Qbase[i,1])*" "*string(Qbase[i,2])*" "*string(Qbase[i,3])*" "*string(Qbase[i,4])*" "*string(T)*"\n")
        end
    end
    println("write "*fff)
end