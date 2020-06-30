using ProgressMeter

function main()
    out_dir="result"
    PARAMDAT="PARAMDAT.json"
    xmax,ymax,nodes,vecAx,vecAy = read_allgrid()
    nt,dt,init_rho,init_u,init_v,init_p,specific_heat_ratio,bdcon,every_outnum,out_file_front,out_ext,restart_file,restart_num = input_para(PARAMDAT)

    Qbase,cellxmax,cellymax,restart_num = set_initQbase(xmax,ymax,restart_file,init_rho,init_u,init_v,init_p,specific_heat_ratio,out_file_front,out_ext,out_dir,restart_num)

    # main loop
    println("threads num ")
    println(Threads.nthreads())
    
    prog = Progress(nt,1)
    for k in 1:nt
        next!(prog)

        volume = set_volume(nodes,cellxmax,cellymax)

        Qbase = set_boundary(Qbase,cellxmax,cellymax,vecAx,vecAy,bdcon)
        Qcon = base_to_conservative(Qbase,cellxmax,cellymax,specific_heat_ratio)
        Qcon_hat = setup_Qcon_hat(Qcon,cellxmax,cellymax,volume)


        cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus = setup_cell_flux_hat(Qbase,Qcon,cellxmax,cellymax,specific_heat_ratio,vecAx,vecAy)

        #println(cell_E_hat_plas[50,:])
        #println(cell_E_hat_plas[51,:])

        E_hat,F_hat = FVS(cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus,Qbase,Qcon,cellxmax,cellymax,specific_heat_ratio,vecAx,vecAy)
        
        #println(E_hat[50,:])
        #println(E_hat[51,:])
        #println(F_hat[:])
        
        RHS = setup_RHS(cellxmax,cellymax,E_hat,F_hat)

        #println(RHS[50,:])
        #println(RHS[51,:])
        
        Qcon_hat = time_integration_explicit(dt,Qcon_hat,RHS,cellxmax,cellymax)
        
        #println(Qcon)
        Qcon = Qhat_to_Q(Qcon_hat,cellxmax,cellymax,volume)

        #println(vecE_hat)
        #println(vecF_hat)
        #println(Qcon)
        
        #throw(UndefVarError(:x))

        Qbase = conservative_to_base(Qcon,cellxmax,cellymax,specific_heat_ratio)

        evalnum=k+restart_num
        if round(evalnum) % every_outnum == 0
            println("nt_______________________________"*string(round(evalnum)))
            output_result(evalnum,Qbase,cellxmax,cellymax,specific_heat_ratio,out_file_front,out_ext,out_dir)
        end


        if Qcon[2,2,4] < 0
            println("pressure are minus !!")
            break
        end

    end
end

# -- main --
main()