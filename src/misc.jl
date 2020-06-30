function set_initQbase(xmax,ymax,restart_file,init_rho,init_u,init_v,init_p,specific_heat_ratio,out_file_front,out_ext,out_dir,restart_num)
    Qbase=[]
    cellxmax = xmax - 1
    cellymax = ymax - 1

    restart_check=0
    try Qbase=setup_restart_value(cellxmax,cellymax,out_dir,restart_file)
        println("Restart "*restart_file)
        restart_check=2
    catch 
        restart_check=1
    end

    if restart_check == 1
        Qbase=setup_init_value(cellxmax,cellymax,init_rho,init_u,init_v,init_p)
        println("Start Initial condition")
        restart_num=0
        output_result(0,Qbase,cellxmax,cellymax,specific_heat_ratio,out_file_front,out_ext,out_dir)
    end

    return Qbase,cellxmax,cellymax,restart_num
end

function setup_init_value(cellxmax,cellymax,init_rho,init_u,init_v,init_p)
    Qbase=zeros(cellxmax, cellymax, 4)
    for i in 1:cellxmax
        for j in 1:cellymax
            Qbase[i,j,1]=init_rho
            Qbase[i,j,2]=init_u
            Qbase[i,j,3]=init_v
            Qbase[i,j,4]=init_p
        end
    end
    return Qbase
end

function setup_restart_value(cellxmax,cellymax,out_dir,restart_file)
    Qbase=zeros(cellxmax,cellymax,4)

    skipnum=1
    fff=[]
    open("result/"*restart_file, "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end
    
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            temp = split(fff[i+skipnum]," ")
            for k in 1:4
                Qbase[i,j,k]=parse(Float64,temp[k]) 
            end
        end
    end
    return Qbase
end

function base_to_conservative(Qbase,cellxmax,cellymax,specific_heat_ratio)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """
    Qcon = zeros(cellxmax,cellymax,4)
    Threads.@threads for i in 1:cellxmax
        for j in 1:cellymax
            Qcon[i,j,1]=Qbase[i,j,1]
            Qcon[i,j,2]=Qbase[i,j,1]*Qbase[i,j,2]
            Qcon[i,j,3]=Qbase[i,j,1]*Qbase[i,j,3]
            Qcon[i,j,4]=Qbase[i,j,4]/(specific_heat_ratio-1)+Qbase[i,j,1]*(Qbase[i,j,2]^2+Qbase[i,j,3]^2)/2
        end
    end
    return Qcon
end

function conservative_to_base(Qcon,cellxmax,cellymax,specific_heat_ratio)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """
    Qbase = zeros(cellxmax,cellymax,4)
    Threads.@threads for i in 1:cellxmax
        for j in 1:cellymax
            Qbase[i,j,1]=Qcon[i,j,1]
            Qbase[i,j,2]=Qcon[i,j,2]/Qcon[i,j,1]
            Qbase[i,j,3]=Qcon[i,j,3]/Qcon[i,j,1]
            Qbase[i,j,4]=(Qcon[i,j,4]-Qcon[i,j,1]*(Qbase[i,j,2]^2+Qbase[i,j,3]^2)/2)*(specific_heat_ratio-1)
        end
    end
    return Qbase
end