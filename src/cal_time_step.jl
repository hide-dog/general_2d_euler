function time_integration_explicit(dt,Qcon_hat,RHS,cellxmax,cellymax)
    Threads.@threads for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 1:4
                Qcon_hat[i,j,k] = Qcon_hat[i,j,k] + dt*RHS[i,j,k]
            end
        end
    end
    return Qcon_hat
end 

function one_wave(Qbase,Qcon,cellxmax,cellymax,vecAx,vecAy,specific_heat_ratio,volume)
    # セル中心のonewave
    A_adv_hat_p = zeros(cellxmax,cellymax,4,4)
    A_adv_hat_m = zeros(cellxmax,cellymax,4,4)
    B_adv_hat_p = zeros(cellxmax,cellymax,4,4)
    B_adv_hat_m = zeros(cellxmax,cellymax,4,4)
    A_beta_sig = zeros(cellxmax,cellymax)
    B_beta_sig = zeros(cellxmax,cellymax)
    beta = 1.1

    for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:2 #A,B
                if k == 1
                    # xi_x = vecAv?
                    kx_av = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1])
                    ky_av = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2])
                elseif k == 2
                    kx_av = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1])
                    ky_av = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2])
                end

                jacob_temp = zeros(4,4)

                rho = Qbase[i,j,1]
                u = Qbase[i,j,2]
                v = Qbase[i,j,3]
                e = Qcon[i,j,4]
                p = Qbase[i,j,4]

                g = specific_heat_ratio
                
                Z = kx_av*u+ky_av*v
                q2   = u^2+v^2
                b1c2 = 0.5*q2*(g-1)
                gebyrho = g*e/rho
                c = (g*rho/p)^0.5

                jacob_temp[1,1] = 0.0
                jacob_temp[1,2] = kx_av
                jacob_temp[1,3] = ky_av
                jacob_temp[1,4] = 0.0
            
                jacob_temp[2,1] = -u*Z + kx_av*b1c2
                jacob_temp[2,2] = Z - (g-2)*kx_av*u
                jacob_temp[2,3] = ky_av*u - kx_av*(g-1)*v
                jacob_temp[2,4] = (g-1)*kx_av
            
                jacob_temp[3,1] = -v*Z + ky_av*b1c2
                jacob_temp[3,2] = kx_av*v - ky_av*(g-1)*u
                jacob_temp[3,3] = Z - ky_av*(g-2)*v
                jacob_temp[3,4] = (g-1)*ky_av
            
                jacob_temp[4,1] = Z*(-gebyrho + 2*b1c2)
                jacob_temp[4,2] = kx_av*(gebyrho-b1c2) - (g-1)*u*Z
                jacob_temp[4,3] = ky_av*(gebyrho-b1c2) - (g-1)*v*Z
                jacob_temp[4,4] = g*Z
                
                sigma = abs(Z) + c*(kx_av^2+ky_av^2)^0.5

                I_temp = zeros(4,4)
                for l in 1:4
                    I_temp[l,l] = beta * sigma
                end

                if k == 1
                    for l in 1:4
                        for m in 1:4
                            A_adv_hat_p[i,j,l,m] = 0.5*(jacob_temp[l,m] + I_temp[l,m])
                            A_adv_hat_m[i,j,l,m] = 0.5*(jacob_temp[l,m] - I_temp[l,m])
                        end
                    end
                    A_beta_sig[i,j] = beta * sigma
                elseif k ==2
                    for l in 1:4
                        for m in 1:4
                            B_adv_hat_p[i,j,l,m] = 0.5*(jacob_temp[l,m] + I_temp[l,m])
                            B_adv_hat_m[i,j,l,m] = 0.5*(jacob_temp[l,m] - I_temp[l,m])
                        end
                    end
                    B_beta_sig[i,j] = beta * sigma
                end
            end
        end
    end
    #=

    println("\n ---------------------------------- \n")
    println(" Aadv ")
    println(A_adv_hat_p[3,3,:,:])
    println("\n ---------------------------------- \n")

    println("\n ---------------------------------- \n")
    println(" A_bs ")
    println(A_beta_sig[3,3])
    println("\n ---------------------------------- \n")
    =#

    return A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_sig, B_beta_sig
end

function set_LDU(dt,Qcon_hat,A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_sig, B_beta_sig,RHS,cellxmax,cellymax,volume)
    D = zeros(cellxmax,cellymax)
    Lx = zeros(cellxmax,cellymax,4,4)
    Ly = zeros(cellxmax,cellymax,4,4)
    Ux = zeros(cellxmax,cellymax,4,4)
    Uy = zeros(cellxmax,cellymax,4,4)

    for i in 1:cellxmax
        for j in 1:cellymax
            D[i,j] = volume[i,j] + dt*(A_beta_sig[i,j] + B_beta_sig[i,j])
            for l in 1:4
                for m in 1:4
                    Lx[i,j,l,m] = -dt*(A_adv_hat_p[i,j,l,m])
                    Ly[i,j,l,m] = -dt*(B_adv_hat_p[i,j,l,m])
                end
            end
            for l in 1:4
                for m in 1:4
                    Ux[i,j,l,m] = dt*(A_adv_hat_m[i,j,l,m])
                    Uy[i,j,l,m] = dt*(B_adv_hat_m[i,j,l,m])
                end
            end
        end
    end

    return Lx,Ly,Ux,Uy,D
end

function lusgs(dt,RHS,cellxmax,cellymax,volume,Lx,Ly,Ux,Uy,D,Delta_Qcon_hat)
    #=
    println("\n ---------------------------------- \n")
    println(" Lx ")
    println(Lx[3,3,:,:])
    println("\n ---------------------------------- \n")
    =#
    
    I = zeros(4,4)
    for i in 1:4
        I[i,i] = 1.0
    end

    for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:4
                RHS[i,j,k] = dt*RHS[i,j,k]
            end
        end
    end   

    # lower sweep
    LdQ = zeros(4)
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for l in 1:4
                for m in 1:4
                    LdQ[l] = Lx[i-1,j,l,m]*Delta_Qcon_hat[i-1,j,m] + Ly[i,j-1,l,m]*Delta_Qcon_hat[i,j-1,m]        
                end
            end

            for l in 1:4
                Delta_Qcon_hat[i,j,l] = D[i,j]^(-1) * (LdQ[l]+RHS[i,j,l])
            end
        end
    end
    #println(Delta_Qcon_hat)

    # upepr sweep
    UdQ = zeros(4)
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for l in 1:4
                for m in 1:4
                    UdQ[l] = Ux[i+1,j,l,m]*Delta_Qcon_hat[i+1,j,m] + Uy[i,j+1,l,m]*Delta_Qcon_hat[i,j+1,m]        
                end
            end

            for l in 1:4
                Delta_Qcon_hat[i,j,l] = Delta_Qcon_hat[i,j,l] - D[i,j]^(-1) * UdQ[l]
            end
        end
    end

    return Delta_Qcon_hat
end

function set_LDU_inner(dt,cell_A_plas,cell_A_minus,cell_B_plas,cell_B_minus,cellxmax,cellymax,volume)
    Dinv  = zeros(cellxmax,cellymax,4,4)
    Lx = zeros(cellxmax,cellymax,4,4)
    Ly = zeros(cellxmax,cellymax,4,4)
    Ux = zeros(cellxmax,cellymax,4,4)
    Uy = zeros(cellxmax,cellymax,4,4)

    I = zeros(4,4)
    for i in 1:4
        I[i,i] = 1.0
    end
    
    D  = zeros(4,4)
    for i in 1:cellxmax
        for j in 1:cellymax
            for l in 1:4
                for m in 1:4
                    Lx[i,j,l,m] = -dt*(cell_A_plas[i,j,l,m])
                    Ly[i,j,l,m] = -dt*(cell_B_plas[i,j,l,m])

                    D[l,m] = volume[i,j]*I[l,m] 
                            + dt*(cell_A_plas[i,j,l,m] + cell_A_minus[i,j,l,m])
                            + dt*(cell_B_plas[i,j,l,m] + cell_B_minus[i,j,l,m])
    
                    Ux[i,j,l,m] = dt*(cell_A_minus[i,j,l,m])
                    Uy[i,j,l,m] = dt*(cell_A_minus[i,j,l,m])
                end
            end

            temp = inv_matrix(D)
            for l in 1:4
                for m in 1:4
                    Dinv[i,j,l,m] = temp[l,m]
                end
            end
        end
    end

    return Lx,Ly,Ux,Uy,Dinv
end

function inv_matrix(D)
    # D = n*n行列
    dim = size(D)[1] 
    # Sweep method
    for l in 1:dim
        aa       = 1.0/D[l,l]
        D[l,l] = 1.0

        for m in 1:dim
            D[l,m] = D[l,m] * aa
        end

        for m in 1:dim
            if m != l
                bb       = D[m,l]
                D[m,l] = 0.0
                for n in 1:dim
                D[m,n] = D[m,n] - bb*D[l,n]
                end
            end
        end 
    end
    return D
end

function inner_iteration(dt,RHS,cellxmax,cellymax,volume,Lx,Ly,Ux,Uy,Dinv,Delta_Qcon_hat)
    
    # lower sweep
    LdQ = zeros(4)
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for l in 1:4
                for m in 1:4
                    LdQ[l] = Lx[i-1,j,l,m]*Delta_Qcon_hat[i-1,j,m] + Ly[i,j-1,l,m]*Delta_Qcon_hat[i,j-1,m]        
                end
            end

            # (LdQ[l]+RHS[i,j,l])
            for l in 1:4
                LdQ[l] = LdQ[l] + RHS[i,j,l]
            end

            # Delta_Qcon_hat[i,j,l] = D[i,j,l,m]^(-1) * (LdQ[m]+RHS[i,j,m])
            for l in 1:4
                temp = 0.0
                for m in 1:4
                    temp = temp + Dinv[i,j,l,m] * LdQ[m]
                end
                Delta_Qcon_hat[i,j,l] = temp
            end
        end
    end
    #println(Delta_Qcon_hat)

    # upepr sweep
    UdQ = zeros(4)
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for l in 1:4
                for m in 1:4
                    UdQ[l] = Ux[i+1,j,l,m]*Delta_Qcon_hat[i+1,j,m] + Uy[i,j+1,l,m]*Delta_Qcon_hat[i,j+1,m]        
                end
            end

            # Delta_Qcon_hat[i,j,l] = Delta_Qcon_hat[i,j,l] - D[i,j,l,m]^(-1) * UdQ[m]
            for l in 1:4
                temp = 0.0
                for m in 1:4
                    temp = temp + Dinv[i,j,l,m] * UdQ[m]
                end
                Delta_Qcon_hat[i,j,l] = Delta_Qcon_hat[i,j,l] - temp
            end
        end
    end 
    return Delta_Qcon_hat
end