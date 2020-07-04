using ProgressMeter

function main()
    out_dir="result"
    PARAMDAT="PARAMDAT.json"
    xmax,ymax,nodes,vecAx,vecAy = read_allgrid()
    nt,dt,init_rho,init_u,init_v,init_p,specific_heat_ratio,bdcon,every_outnum,out_file_front,out_ext,restart_file,restartnum,time_integ,init_small,norm_ok,nt_lusgs = input_para(PARAMDAT)

    Qbase,cellxmax,cellymax,restartnum = set_initQbase(xmax,ymax,restart_file,init_rho,init_u,init_v,init_p,specific_heat_ratio,out_file_front,out_ext,out_dir,restartnum)

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

        #println(Qbase[:,1,2])
        #println(Qbase[:,1,3])

        cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus = setup_cell_flux_hat(Qbase,Qcon,cellxmax,cellymax,specific_heat_ratio,vecAx,vecAy)

        E_hat,F_hat = FVS(cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus,Qbase,Qcon,cellxmax,cellymax,specific_heat_ratio,vecAx,vecAy)
        
        RHS = setup_RHS(cellxmax,cellymax,E_hat,F_hat)
        
        # implicit scheme
        #= A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_shig, B_beta_shig = one_wave(Qbase,Qcon,cellxmax,cellymax,vecAx,vecAy,specific_heat_ratio,volume)
        nt_lusgs = Int(1.0*10^(5))
        for i in 1:nt_lusgs
            Qcon_hat,Lx,Ly,Ux,Uy = lusgs(dt,Qcon_hat,A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_shig, B_beta_shig,jalphaP, jbetaP,RHS,cellxmax,cellymax,volume)
            res = set_res(Qcon_hat,Lx,Ly,Ux,Uy,RHS,cellxmax,cellymax)
            norm2 = check_converge(res,RHS,cellxmax,cellymax)
            if norm2 < norm_ok
                break
            end 
            if i == nt_lusgs
                println("\n ---------------------------------- \n")
                println(" i == nt_lusgs ! \n")
                println(" check lusgs")
                println("\n ---------------------------------- \n")
            end 
            
        end
        =#
        
        evalnum = k+restartnum
        
        if time_integ == "1"
            # exlicit scheme
            Qcon_hat = time_integration_explicit(dt,Qcon_hat,RHS,cellxmax,cellymax)
        elseif time_integ == "2"
            # implicit scheme
            A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_shig, B_beta_shig = one_wave(Qbase,Qcon,cellxmax,cellymax,vecAx,vecAy,specific_heat_ratio,volume)
            Qcon_hat,Lx,Ly,Ux,Uy = lusgs(dt,Qcon_hat,A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_shig, B_beta_shig,RHS,cellxmax,cellymax,volume)
            res = set_res(Qcon_hat,Lx,Ly,Ux,Uy,RHS,cellxmax,cellymax)
            norm2 = check_converge(res,RHS,cellxmax,cellymax,init_small)
            println("\n ---------------------------------- \n")
            println("nt : "*string(round(evalnum)))
            println("density res:"*string(norm2[1]) * "  energy res:"*string(norm2[4]))
            println("\n ---------------------------------- \n")
        else
            println(time_integ)
            println("pressure check time_integ")
            throw(UndefVarError(:x))
        end
        
        Qcon = Qhat_to_Q(Qcon_hat,cellxmax,cellymax,volume)
        Qbase = conservative_to_base(Qcon,cellxmax,cellymax,specific_heat_ratio)

        
        if round(evalnum) % every_outnum == 0
            println("nt_______________________________"*string(round(evalnum)))
            output_result(evalnum,Qbase,cellxmax,cellymax,specific_heat_ratio,out_file_front,out_ext,out_dir)
        end


        if Qcon[2,2,4] < 0
            println("pressure are minus !!")
            break
        end
        if Qcon[2,2,1] == NaN
            println(" NaN !!")
            break
        end

    end
end

# -- main --
main()