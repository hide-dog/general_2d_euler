using ProgressMeter

function main()
    out_dir="result"
    PARAMDAT="PARAMDAT.json"
    xmax,ymax,nodes,vecAx,vecAy = read_allgrid()
    nt,dt,init_rho,init_u,init_v,init_p,specific_heat_ratio,bdcon,every_outnum,out_file_front,out_ext,restart_file,restartnum,time_integ,init_small,norm_ok,nt_lusgs = input_para(PARAMDAT)

    Qbase,cellxmax,cellymax,restartnum = set_initQbase(xmax,ymax,restart_file,init_rho,init_u,init_v,init_p,specific_heat_ratio,out_file_front,out_ext,out_dir,restartnum)

    # init Delta_Qcon_hat
    Delta_Qcon_hat = zeros(cellxmax,cellymax,4)
    
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
        
        cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus = cal_jacobi(Qbase,Qcon,cellxmax,cellymax,specific_heat_ratio,vecAx,vecAy)
        cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus = setup_cell_flux_hat(Qcon,cellxmax,cellymax,cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus)

        println(cell_E_hat_plas[2,2,:])
        println(cell_E_hat_minus[2,2,:])

        E_hat,F_hat = FVS(cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus,cellxmax,cellymax)
        
        RHS = setup_RHS(cellxmax,cellymax,E_hat,F_hat,Qcon_hat)
        evalnum = k+restartnum

        if time_integ == "1"
            # exlicit scheme
            Qcon_hat = time_integration_explicit(dt,Qcon_hat,RHS,cellxmax,cellymax)
        elseif time_integ == "2"
            # implicit scheme
            A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_sig, B_beta_sig = one_wave(Qbase,Qcon,cellxmax,cellymax,vecAx,vecAy,specific_heat_ratio,volume)
            Lx,Ly,Ux,Uy,D = set_LDU(dt,Qcon_hat,A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_sig, B_beta_sig,RHS,cellxmax,cellymax,volume)

            Delta_Qcon_hat_temp = zeros(cellxmax,cellymax,4)
            for lusgs_t in 1:nt_lusgs
                Delta_Qcon_hat = lusgs(dt,RHS,cellxmax,cellymax,volume,Lx,Ly,Ux,Uy,D,Delta_Qcon_hat)
                
                res = set_res(Delta_Qcon_hat,Delta_Qcon_hat_temp,cellxmax,cellymax)
                norm2 = check_converge(res,RHS,cellxmax,cellymax,init_small)
                
                Delta_Qcon_hat_temp = copy(Delta_Qcon_hat)
                
                if norm2[1] < norm_ok && norm2[4] < norm_ok
                    println("\n ---------------------------------- \n")
                    println("nt : "*string(round(evalnum)))
                    println("density res:"*string(norm2[1]) * "  energy res:"*string(norm2[4]))
                    println("\n ---------------------------------- \n")
                    break
                end 
                
            end

            # lusgs = update
            for i in 2:cellxmax-1
                for j in 2:cellymax-1
                    for l in 1:4
                        Qcon_hat[i,j,l] = Qcon_hat[i,j,l] + Delta_Qcon_hat[i,j,l]
                    end
                end
            end

            # 内部反復
            #=
            nt_inner = 10
            Qcon_hat_n = copy(Qcon_hat)
            Qcon_hat_m1 = copy(Qcon_hat)
            Delta_Qcon_hat_m = copy(Delta_Qcon_hat)

            for inner_t in 1:nt_inner
                Delta_Qcon_hat_temp = zeros(cellxmax,cellymax,4)

                for lusgs_t in 1:nt_lusgs
                    Qcon_m1 = Qhat_to_Q(Qcon_hat_m1,cellxmax,cellymax,volume)
                    Qbase_m1 = conservative_to_base(Qcon_m1,cellxmax,cellymax,specific_heat_ratio)

                    cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus = cal_jacobi(Qbase_m1,Qcon_m1,cellxmax,cellymax,specific_heat_ratio,vecAx,vecAy)
                    cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus = setup_cell_flux_hat(Qcon_m1,cellxmax,cellymax,cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus)
                    E_hat,F_hat = FVS(cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus,cellxmax,cellymax)
                    
                    # lusgsと違うので注意
                    RHS = set_RHS_inner(cellxmax,cellymax,dt,Qcon_hat_n,Qcon_hat_m1,E_hat,F_hat)
                    Lx,Ly,Ux,Uy,Dinv = set_LDU_inner(dt,cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus,cellxmax,cellymax,volume)

                    Delta_Qcon_hat_m = inner_iteration(dt,RHS,cellxmax,cellymax,volume,Lx,Ly,Ux,Uy,Dinv,Delta_Qcon_hat_m)


                    res = set_res(Delta_Qcon_hat,Delta_Qcon_hat_temp,cellxmax,cellymax)
                    norm2 = check_converge(res,RHS,cellxmax,cellymax,init_small)
                    Delta_Qcon_hat_temp = copy(Delta_Qcon_hat_m)

                    if norm2[1] < norm_ok && norm2[4] < norm_ok
                        println("\n ---------------------------------- \n")
                        println("nt : "*string(round(inner_t)))
                        println("density res:"*string(norm2[1]) * "  energy res:"*string(norm2[4]))
                        println("\n ---------------------------------- \n")
                        break
                    end
                end

                # inner update
                for i in 2:cellxmax-1
                    for j in 2:cellymax-1
                        for l in 1:4
                            Qcon_hat_m1[i,j,l] = Qcon_hat_m1[i,j,l] + Delta_Qcon_hat_m[i,j,l]
                        end
                    end
                end
            end

            Qcon_hat = copy(Qcon_hat_m1)
            =#
        
            #println(Qcon_hat[:,:,1])
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