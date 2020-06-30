using ProgressMeter

function main()
    out_dir="result"
    PARAMDAT="PARAMDAT.json"
    nodes,elements,ele_bdx,ele_bdy,vecAx,vecAy,num_cell,cell_Enum,cell_Fnum,volume,fx_volume,fy_volume,cell_vecAx,cell_vecAy=read_allgrid()
    nt,dt,init_rho,init_u,init_v,init_p,specific_heat_ratio,bdcon,every_outnum,out_file_front,out_ext,restart_file,restart_num=input_para(PARAMDAT)

    Qbase=[]
    Qcon=[]

    restart_check=0
    try Qbase=setup_restart_value(num_cell,out_dir,restart_file)
        println("Restart "*restart_file)
        restart_check=2
    catch 
        restart_check=1
    end

    if restart_check == 0
        Qbase=setup_init_value(num_cell,init_rho,init_u,init_v,init_p)
        println("Start Initial condition")
        restart_num=0
        output_result(0,Qbase,specific_heat_ratio,out_file_front,out_ext,out_dir)
    elseif restart_check == 1
        println("restart error")
        println("Please check restart file")
        throw(UndefVarError(:x))
    end 



    # main loop
    prog = Progress(nt,1)
    for k in 1:nt
        next!(prog)
        Qcon=base_to_conservative(Qbase,specific_heat_ratio)

        cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus=setup_cell_flux_hat(Qbase,Qcon,specific_heat_ratio,cell_vecAx,cell_vecAy)

        #println(cell_E_hat_plas[50,:])
        #println(cell_E_hat_plas[51,:])

        Qcon_hat=setup_Qcon_hat(Qcon,volume)

        E_hat,F_hat=FVS(cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus,Qbase,Qcon,vecAx,vecAy,ele_bdx,ele_bdy,specific_heat_ratio,bdcon,cell_Enum,cell_Fnum)
        
        #println(E_hat[50,:])
        #println(E_hat[51,:])
        #println(F_hat[:])
        
        RHS=setup_RHS(cell_Enum,cell_Fnum,E_hat,F_hat)

        #println(RHS[50,:])
        #println(RHS[51,:])
        
        Qcon_hat=time_integration_explicit(dt,Qcon_hat,RHS)
        
        #println(Qcon)
        Qcon=Qhat_to_Q(Qcon_hat,volume)

        #println(vecE_hat)
        #println(vecF_hat)
        #println(Qcon)
        
        #throw(UndefVarError(:x))

        Qbase=conservative_to_base(Qcon,specific_heat_ratio)

        evalnum=k+restart_num
        if round(evalnum) % every_outnum == 0
            println("nt_______________________________"*string(round(evalnum)))
            output_result(evalnum,Qbase,specific_heat_ratio,out_file_front,out_ext,out_dir)
        end


        if Qcon[size(Qcon)[1],size(Qcon)[2]] < 0
            println("pressure are minus !!")
            break
        end

    end
end

# -- main --
main()