function set_volume(nodes,cellxmax,cellymax)
    volume = zeros(cellxmax,cellymax)
    for i in 1:cellxmax
        for j in 1:cellymax
            vec_r1x = nodes[i+1,j+1,1] - nodes[i,j,1]
            vec_r1y = nodes[i+1,j+1,2] - nodes[i,j,2]
            vec_r2x = nodes[i,j+1,1] - nodes[i+1,j,1]
            vec_r2y = nodes[i,j+1,2] - nodes[i+1,j,2]

            volume[i,j] = abs(vec_r1x*vec_r2y - vec_r1y*vec_r2x) /2
        end
    end
    return volume
end 

function setup_cell_flux_hat(Qcon,cellxmax,cellymax,cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus)    
    cell_E_hat_plas  = zeros(cellxmax,cellymax,4)
    cell_E_hat_minus = zeros(cellxmax,cellymax,4)
    cell_F_hat_plas  = zeros(cellxmax,cellymax,4)
    cell_F_hat_minus = zeros(cellxmax,cellymax,4)

    for i in 1:cellxmax
        for j in 1:cellymax
            for l in 1:4
                for m in 1:4
                    # Ehat = Ahat * Qhat = Ahat * Q/J = 1/J*Ahat * Q
                    # 1/J*Ahatの計算をvecAを使ったため，ここではQconを掛けてる
                    cell_E_hat_plas[i,j,l]  += cell_Ahat_plas[i,j,l,m]*Qcon[i,j,m]
                    cell_E_hat_minus[i,j,l] += cell_Ahat_minus[i,j,l,m]*Qcon[i,j,m]
                    cell_F_hat_plas[i,j,l]  += cell_Bhat_plas[i,j,l,m]*Qcon[i,j,m]
                    cell_F_hat_minus[i,j,l] += cell_Bhat_minus[i,j,l,m]*Qcon[i,j,m]
                end
            end
        end
    end

    return cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus
end

function cal_jacobi(Qbase,Qcon,cellxmax,cellymax,specific_heat_ratio,vecAx,vecAy)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    
    q=u^2+v^2
    c=(rp/rho)^0.5
    theta=Z=U= A^xi_x * u   +   A^xi_y * v
    H=(e+p)/rho
    b1=q^2/2 * (r-1)/(c^2)
    b2=(r-1)/(c^2)

    xi_x_bar= xi_x / (xi_x^2+xi_y^2)^0.5
    xi_y_bar= xi_y / (xi_x^2+xi_y^2)^0.5
    U_bar=Z_bar=theta_bar=


    A = Rleft * Lambda * Rright
    """

    cell_Ahat_plas  = zeros(cellxmax,cellymax,4,4)
    cell_Ahat_minus = zeros(cellxmax,cellymax,4,4)
    cell_Bhat_plas  = zeros(cellxmax,cellymax,4,4)
    cell_Bhat_minus = zeros(cellxmax,cellymax,4,4)
    r=specific_heat_ratio

    for i in 1:cellxmax
        for j in 1:cellymax
            q=(Qbase[i,j,2]^2+Qbase[i,j,3]^2)^0.5
            c=0
            try c=(r*Qbase[i,j,4]/Qbase[i,j,1])^0.5
            catch
                println("\n"*string(i)*","*string(j)*" Qbase error")
                println(Qbase[i,j,:])
                println("\n")
                throw(UndefVarError(:x))
            end

            H  = (Qcon[i,j,4]+Qbase[i,j,4])/Qbase[i,j,1]
            b1 = q^2/2 * (r-1)/(c^2)
            b2 = (r-1)/(c^2)
            u  = Qbase[i,j,2]
            v  = Qbase[i,j,3]

            k_x = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1])
            k_y = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2])
            
            Rleft_E,Lambda_E,Rright_E = Eigenvalue_vector(q,c,H,b1,b2,u,v,k_x,k_y)
            A_plas, A_minus = RLmbdaR(Rleft_E,Lambda_E,Rright_E)

            for l in 1:4
                for m in 1:4
                    cell_Ahat_plas[i,j,l,m]  = A_plas[l,m]
                    cell_Ahat_minus[i,j,l,m] = A_minus[l,m]
                end
            end

            k_x = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1])
            k_y = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2])

            Rleft_E,Lambda_E,Rright_E = Eigenvalue_vector(q,c,H,b1,b2,u,v,k_x,k_y)
            A_plas, A_minus = RLmbdaR(Rleft_E,Lambda_E,Rright_E)

            for l in 1:4
                for m in 1:4
                    cell_Bhat_plas[i,j,l,m]  = A_plas[l,m]
                    cell_Bhat_minus[i,j,l,m] = A_minus[l,m]
                end
            end
        end
    end

    return cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus
end


function Eigenvalue_vector(q,c,H,b1,b2,u,v,k_x,k_y)
    Z= k_x*u + k_y*v
    k_x_bar=k_x / (k_x^2+k_y^2)^0.5
    k_y_bar=k_y / (k_x^2+k_y^2)^0.5
    Z_bar=Z / (k_x^2+k_y^2)^0.5
    
    Lambda = [Z-c*(k_x^2+k_y^2)^0.5 
              Z
              Z+c*(k_x^2+k_y^2)^0.5 
              Z]

    Rleft = [[1 1 1 0]
             [u-k_x_bar*c u u+k_x_bar*c -k_y_bar]
             [v-k_y_bar*c v v+k_y_bar*c k_x_bar]
             [H-c*Z_bar q^2/2 H+c*Z_bar -(k_y_bar*u-k_x_bar*v)]]

    Rright = [[(b1+Z_bar/c)/2 -(k_x_bar/c+b2*u)/2 -(k_y_bar/c+b2*v)/2 b2/2]
              [1-b1 b2*u b2*v -b2]
              [(b1-Z_bar/c)/2 (k_x_bar/c-b2*u)/2 (k_y_bar/c-b2*v)/2 b2/2]
              [k_y_bar*u-k_x_bar*v -k_y_bar k_x_bar 0]]
    
    return Rleft,Lambda,Rright
end

function RLmbdaR(R,Lam,Rm)
    Lamp=zeros(4,4)
    Lamm=zeros(4,4)
    for i in 1:length(Lam)
        Lamp[i,i]=(Lam[i]+abs(Lam[i]))/2
        Lamm[i,i]=(Lam[i]-abs(Lam[i]))/2
    end
    Ap=nn_inner_product(nn_inner_product(R,Lamp),Rm)
    Am=nn_inner_product(nn_inner_product(R,Lamm),Rm)
    return Ap,Am
end

function nn_inner_product(a,b)
    temp=zeros(size(a)[1],size(a)[1])
    for i in 1:size(a)[1]
        for j in 1:size(a)[1]
            for k in 1:size(a)[1]
                temp[i,j] += a[i,k]*b[k,j]
            end
        end
    end
    return temp    
end


function FVS(cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus,cellxmax,cellymax)
    E_hat = zeros(cellxmax+1,cellymax,4)
    F_hat = zeros(cellxmax,cellymax+1,4)

    Threads.@threads for i in 2:cellxmax+1-1
        for j in 2:cellymax-1
            for k in 1:4
                E_hat[i,j,k] = cell_E_hat_plas[i-1,j,k] + cell_E_hat_minus[i,j,k]
            end
        end
    end
    Threads.@threads for i in 2:cellxmax-1
        for j in 2:cellymax+1-1
            for k in 1:4
                F_hat[i,j,k] = cell_F_hat_plas[i,j-1,k] + cell_F_hat_minus[i,j,k]
            end
        end
    end

    return E_hat,F_hat
end

function setup_RHS(cellxmax,cellymax,E_hat,F_hat,Qcon_hat)
    RHS = zeros(cellxmax,cellymax,4)
    
    Threads.@threads for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 1:4
                RHS[i,j,k] = -(E_hat[i+1,j,k]-E_hat[i,j,k]+F_hat[i,j+1,k]-F_hat[i,j,k])
            end
        end
    end
    return RHS
end

function set_RHS_inner(cellxmax,cellymax,dt,Qcon_hat_n,Qcon_hat_m1,E_hat,F_hat)
    RHS = zeros(cellxmax,cellymax,4)
    
    Threads.@threads for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 1:4
                RHS[i,j,k] = -(Qcon_hat_m1[i,j,k]-Qcon_hat_n[i,j,k]) -dt*(E_hat[i+1,j,k]-E_hat[i,j,k]+F_hat[i,j+1,k]-F_hat[i,j,k])
            end
        end
    end
    return RHS
end

function setup_Qcon_hat(Qcon,cellxmax,cellymax,volume)
    Qcon_hat=ones(cellxmax,cellymax,4)
    Threads.@threads for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:4
                Qcon_hat[i,j,k] = Qcon[i,j,k]*volume[i,j]
            end
        end
    end
    return Qcon_hat
end

function Qhat_to_Q(Qcon_hat,cellxmax,cellymax,volume)
    Qcon = ones(cellxmax,cellymax,4)
    Threads.@threads for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:4
                Qcon[i,j,k] = Qcon_hat[i,j,k]/volume[i,j]
            end
        end
    end
    return Qcon
end