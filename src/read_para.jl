import JSON

function make_json(PARAMDAT)
    fff=[]
    open(PARAMDAT, "r") do f
        fff=read(f,String)
    end

    fff=replace(fff,r"#(.+)\n" => "\n")
    
    read_PARAMDAT="read_"*PARAMDAT
    open(read_PARAMDAT, "w") do f
        write(f,fff)
    end
    return read_PARAMDAT
end 

function read_json(read_PARAMDAT)
    dict=1
    open(read_PARAMDAT, "r") do f
        dicttxt = read(f,String)  # file information to string
        dict=JSON.parse(dicttxt)  # parse and transform data
    end
    return dict
end

function read_para(dict)
    restartnum=dict["restartnum"]
    Restart_file=dict["Restart"]*restartnum*dict["in_ext"]
    restartnum=Int(parse(Float64,restartnum))

    nt=Int(parse(Float64,dict["nt"]))                    # 時間ステップ数
    dt=parse(Float64,dict["dt"])                    # 時間刻み幅
    
    # 初期値
    init_rho=parse(Float64,dict["init_rho"])
    init_u=parse(Float64,dict["init_u"])
    init_v=parse(Float64,dict["init_v"])
    init_p=parse(Float64,dict["init_p"])

    specific_heat_ratio=parse(Float64,dict["specific_heat_ratio"])

    # 境界条件
    #= 
    bdcon[i][j]
    i:境界条件番号(1,2,3..)
    j=1:境界条件番号（0or1）
    j=2:rho
    j=3:u
    j=4:v
    j=5:p
    =#
    bdcon=[]
    k=1
    while true
        try
            bd=[parse(Int,dict["bd"*string(k)*"_con"]),parse(Float64,dict["bd"*string(k)*"_rho"]),
                parse(Float64,dict["bd"*string(k)*"_u"]),parse(Float64,dict["bd"*string(k)*"_v"]),
                parse(Float64,dict["bd"*string(k)*"_p"]),parse(Float64,dict["bd"*string(k)*"_T"])]
            push!(bdcon,bd)
        catch
            break
        end
        k+=1
    end
    
    # 出力情報
    every_outnum=Int(parse(Float64,dict["every_outnum"]))
    out_file_front=dict["out_file_front"]
    out_ext=dict["out_ext"]

    return nt,dt,init_rho,init_u,init_v,init_p,specific_heat_ratio,bdcon,every_outnum,out_file_front,out_ext,Restart_file,restartnum
end

function input_para(PARAMDAT)
   read_PARAMDAT=make_json(PARAMDAT)
   dict=read_json(read_PARAMDAT)
   nt,dt,init_rho,init_u,init_v,init_p,specific_heat_ratio,bdcon,every_outnum,out_file_front,out_ext,Restart_file,restartnum=read_para(dict)
   println("fin read para")
   return nt,dt,init_rho,init_u,init_v,init_p,specific_heat_ratio,bdcon,every_outnum,out_file_front,out_ext,Restart_file,restartnum
end