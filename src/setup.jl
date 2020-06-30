function setup_init_value(num_cell,init_rho,init_u,init_v,init_p)
    Qbase=ones(num_cell,4)
    for i in 1:size(Qbase)[1]
        Qbase[i,1]=init_rho
        Qbase[i,2]=init_u
        Qbase[i,3]=init_v
        Qbase[i,4]=init_p
    end
    return Qbase
end

function setup_restart_value(num_cell,indir,restart)
    skipnum=1
    fff=[]
    open(indir*"/"*restart, "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end
    Qbase=ones(num_cell,4)
    
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)-1
            Qbase[i,j]=parse(Float64,temp[j+1]) 
        end
    end
    return Qbase
end

function base_to_conservative(Qbase,specific_heat_ratio)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """
    Qcon=ones(size(Qbase)[1],4)
    for i in 1:size(Qcon)[1]
        Qcon[i,1]=Qbase[i,1]
        Qcon[i,2]=Qbase[i,1]*Qbase[i,2]
        Qcon[i,3]=Qbase[i,1]*Qbase[i,3]
        Qcon[i,4]=Qbase[i,4]/(specific_heat_ratio-1)+Qbase[i,1]*(Qbase[i,2]^2+Qbase[i,3]^2)/2
    end
    return Qcon
end

function conservative_to_base(Qcon,specific_heat_ratio)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """
    Qbase=ones(size(Qcon)[1],4)
    for i in 1:size(Qbase)[1]
        Qbase[i,1]=Qcon[i,1]
        Qbase[i,2]=Qcon[i,2]/Qcon[i,1]
        Qbase[i,3]=Qcon[i,3]/Qcon[i,1]
        Qbase[i,4]=(Qcon[i,4]-Qcon[i,1]*(Qbase[i,2]^2+Qbase[i,3]^2)/2)*(specific_heat_ratio-1)
    end
    return Qbase
end

function Q_to_flux(Qcon,Qbase)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """
    E=ones(size(Qcon)[1],4)
    for i in 1:size(E)[1]
        E[i,1]=Qcon[i,2]
        E[i,2]=Qcon[i,2]*Qbase[i,2]+Qbase[i,4]
        E[i,3]=Qcon[i,2]*Qbase[i,3]
        E[i,4]=(Qcon[i,4]+Qbase[i,4])*Qbase[i,2]
    end

    F=ones(size(Qcon)[1],4)
    for i in 1:size(F)[1]
        F[i,1]=Qcon[i,3]
        F[i,2]=Qcon[i,3]*Qbase[i,2]
        F[i,3]=Qcon[i,3]*Qbase[i,3]+Qbase[i,4]
        F[i,4]=(Qcon[i,4]+Qbase[i,4])*Qbase[i,3]
    end
    return E,F
end


"""
境界における値を定めて，境界におけるフラックスを定義する？
これで境界条件とするので，FVSのときに境界条件を定義
"""

function setup_cell_flux_hat(Qbase,Qcon,specific_heat_ratio,cell_vecAx,cell_vecAy)
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
    cell_E_hat_plas=zeros(size(Qcon)[1],4)
    cell_E_hat_minus=zeros(size(Qcon)[1],4)
    cell_F_hat_plas=zeros(size(Qcon)[1],4)
    cell_F_hat_minus=zeros(size(Qcon)[1],4)

    r=specific_heat_ratio

    for i in 1:size(Qcon)[1]
        q=(Qbase[i,2]^2+Qbase[i,3]^2)^0.5
        c=0
        try c=(r*Qbase[i,4]/Qbase[i,1])^0.5
        catch
            println("\n"*string(i)*"Qbase pressure=minus")
            println(Qbase[i,:])
            throw(UndefVarError(:x))
        end

        H=(Qcon[i,4]+Qbase[i,4])/Qbase[i,1]
        b1=q^2/2 * (r-1)/(c^2)
        b2=(r-1)/(c^2)
        u=Qbase[i,2]
        v=Qbase[i,3]

        k_x=cell_vecAx[i,1]
        k_y=cell_vecAx[i,2]
        
        
        Rleft_E,Lambda_E,Rright_E = Eigenvalue_vector(q,c,H,b1,b2,u,v,k_x,k_y)

        A_plas, A_minus = RLmbdaR(Rleft_E,Lambda_E,Rright_E)

        for j in 1:4
            for k in 1:4
                cell_E_hat_plas[i,j] += A_plas[j,k]*Qcon[i,k]
                cell_E_hat_minus[i,j] += A_minus[j,k]*Qcon[i,k]
            end
        end
"""
        if i==51
            println(Qcon[i-1,:])
            println(Qcon[i,:])

            println(A_plas)
            println(A_minus)

            println(cell_E_hat_plas[i-1,:])
            println(cell_E_hat_minus[i,:])
            throw(UndefVarError(:x))
        end"""

        k_x=cell_vecAy[i,1]
        k_y=cell_vecAy[i,2]
        Rleft_E,Lambda_E,Rright_E = Eigenvalue_vector(q,c,H,b1,b2,u,v,k_x,k_y)
        A_plas, A_minus = RLmbdaR(Rleft_E,Lambda_E,Rright_E)

        for j in 1:4
            for k in 1:4
                cell_F_hat_plas[i,j] += A_plas[j,k]*Qcon[i,k]
                cell_F_hat_minus[i,j] += A_minus[j,k]*Qcon[i,k]
            end
        end  
    end
    return cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus
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


function FVS(E_plas,E_minus,F_plas,F_minus,Qbase,Qcon,vecAx,vecAy,ele_bdx,ele_bdy,specific_heat_ratio,bdcon,cell_Enum,cell_Fnum)
    E_hat=zeros(size(ele_bdx)[1],size(Qbase)[2])

    for i in 1:size(ele_bdx)[1]
        temp=zeros(size(Qbase)[2])
        if ele_bdx[i,4] < 0
            """
            ele_bdx[4]:境界を挟むセル番号
            ele_bdx[5]:境界を挟むセル番号

            wallnum:x=1
                    y=2
            """
            cellnum1=abs(ele_bdx[i,4])
            cellnum2=ele_bdx[i,5]
            bd_E_hat_plas,bd_E_hat_minus,bd_F_hat_plas,bd_F_hat_minus = boundary_condition(bdcon[cellnum1][1],bdcon[cellnum1][2:length(bdcon[cellnum1])],Qbase[cellnum2,:],vecAx[i,3:5],specific_heat_ratio,1)
            temp=bd_E_hat_plas[:]
            for j in 1:4
                temp[j] = temp[j]+E_minus[ele_bdx[i,5],j]
            end
        elseif ele_bdx[i,5] < 0
            cellnum1=abs(ele_bdx[i,5])
            cellnum2=ele_bdx[i,4]
            bd_E_hat_plas,bd_E_hat_minus,bd_F_hat_plas,bd_F_hat_minus = boundary_condition(bdcon[cellnum1][1],bdcon[cellnum1][2:length(bdcon[cellnum1])],Qbase[cellnum2,:],vecAx[i,3:5],specific_heat_ratio,1)
            temp=bd_E_hat_minus[:]
            for j in 1:4
                temp[j] = E_plas[ele_bdx[i,4],j]+temp[j]
            end
        else
            for j in 1:4
                temp[j] = E_plas[ele_bdx[i,4],j]+E_minus[ele_bdx[i,5],j]
            end
        end
        for j in 1:4
            E_hat[i,j]=temp[j]
        end
    end

    F_hat=zeros(size(ele_bdy)[1],size(Qbase)[2])
    for i in 1:size(ele_bdy)[1]
        temp=zeros(size(Qbase)[2])
        if ele_bdy[i,4] < 0
            cellnum1=abs(ele_bdy[i,4])
            cellnum2=ele_bdy[i,5]
            bd_E_hat_plas,bd_E_hat_minus,bd_F_hat_plas,bd_F_hat_minus = boundary_condition(bdcon[cellnum1][1],bdcon[cellnum1][2:length(bdcon[cellnum1])],Qbase[cellnum2,:],vecAy[i,3:5],specific_heat_ratio,2)
            temp=bd_F_hat_plas[:]
            for j in 1:4
                temp[j] = temp[j]+F_minus[ele_bdy[i,5],j]
            end
        elseif ele_bdy[i,5] < 0
            cellnum1=abs(ele_bdy[i,5])
            cellnum2=ele_bdy[i,4]
            bd_E_hat_plas,bd_E_hat_minus,bd_F_hat_plas,bd_F_hat_minus = boundary_condition(bdcon[cellnum1][1],bdcon[cellnum1][2:length(bdcon[cellnum1])],Qbase[cellnum2,:],vecAy[i,3:5],specific_heat_ratio,2)
            temp=bd_F_hat_minus[:]
            for j in 1:4
                temp[j] = F_plas[ele_bdy[i,4],j]+temp[j]
            end
        else
            for j in 1:4
                temp[j] = F_plas[ele_bdy[i,4],j]+F_minus[ele_bdy[i,5],j]
            end
        end
        for j in 1:4
            F_hat[i,j]=temp[j]
        end
    end

    return E_hat,F_hat
end

function boundary_condition(bd_con,bd_val,Qbase,bd_vecA,specific_heat_ratio,wallnum)
    """
    bd_con = bd_num
    bd_val =[rho,u,v,p]
    Qbase =[rho,u,v,p]
    bd_vecA = 境界のvecA vecAx or vecAy

    wallnum:Eを求めるときは2
            Fを求めるときは3
    """
    flux=[]
    
    bd_con=Int(bd_con)
    if bd_con == 0
        bdQbase=bd_val[:]
    elseif bd_con == 1
        bdQbase=Qbase[:]
    elseif bd_con == 2
        bdQbase=wall_Qbase(Qbase,bd_vecA,bd_val,wallnum)
    else
        "change boundary condition"
    end
    bdQcon=bd_base_to_con(bdQbase,specific_heat_ratio)
    bd_E_hat_plas,bd_E_hat_minus,bd_F_hat_plas,bd_F_hat_minus=bd_flux(bdQbase,bdQcon,bd_vecA,specific_heat_ratio)
    #println(bd_E_hat_plas)

    return bd_E_hat_plas,bd_E_hat_minus,bd_F_hat_plas,bd_F_hat_minus
end

function wall_Qbase(Qbase,bd_vecA,bd_val,wallnum)
    """bd_val
    1:rho
    2:U
    3:v
    4:p
    5:T
    """
    bdQbase=zeros(length(Qbase))
    a=bd_vecA[1]
    b=bd_vecA[2]
    u=Qbase[2]
    v=Qbase[3]

    bdQbase[1]=Qbase[1]
    bdQbase[2]=((-a^2+b^2)*u-2*a*b*v)/(a^2+b^2)
    bdQbase[3]=(-2*a*b*u+(a^2-b^2)*v)/(a^2+b^2)
    bdQbase[4]=Qbase[4]  
    
    return bdQbase
end

function bd_base_to_con(bdQbase,specific_heat_ratio)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """
    Qcon=zeros(length(bdQbase))
    
    Qcon[1]=bdQbase[1]
    Qcon[2]=bdQbase[1]*bdQbase[2]
    Qcon[3]=bdQbase[1]*bdQbase[3]
    Qcon[4]=bdQbase[4]/(specific_heat_ratio-1)+bdQbase[1]*(bdQbase[2]^2+bdQbase[3]^2)/2
    
    return Qcon
end

function bd_flux(bdQbase,bdQcon,vecA,specific_heat_ratio)

    bd_E_hat_plas=zeros(length(bdQcon))
    bd_E_hat_minus=zeros(length(bdQcon))
    bd_F_hat_plas=zeros(length(bdQcon))
    bd_F_hat_minus=zeros(length(bdQcon))

    r=specific_heat_ratio

    q=(bdQbase[2]^2+bdQbase[3]^2)^0.5
    c=0
    try c=(r*bdQbase[4]/bdQbase[1])^0.5
    catch
        println("\n bdQbase pressure=minus")
        println(bdQbase[:])
        throw(UndefVarError(:x))
    end

    H=(bdQcon[4]+bdQbase[4])/bdQbase[1]
    b1=q^2/2 * (r-1)/(c^2)
    b2=(r-1)/(c^2)
    u=bdQbase[2]
    v=bdQbase[3]

    k_x=vecA[1]
    k_y=vecA[2]
    Rleft_E,Lambda_E,Rright_E = Eigenvalue_vector(q,c,H,b1,b2,u,v,k_x,k_y)
    A_plas, A_minus = RLmbdaR(Rleft_E,Lambda_E,Rright_E)

    for j in 1:4
        for k in 1:4
            bd_E_hat_plas[j] += A_plas[j,k]*bdQcon[k]
            bd_E_hat_minus[j] += A_minus[j,k]*bdQcon[k]
        end
    end

    k_x=vecA[1]
    k_y=vecA[2]
    Rleft_E,Lambda_E,Rright_E = Eigenvalue_vector(q,c,H,b1,b2,u,v,k_x,k_y)
    A_plas, A_minus = RLmbdaR(Rleft_E,Lambda_E,Rright_E)

    for j in 1:4
        for k in 1:4
            bd_F_hat_plas[j] += A_plas[j,k]*bdQcon[k]
            bd_F_hat_minus[j] += A_minus[j,k]*bdQcon[k]
        end
    end  

    return bd_E_hat_plas,bd_E_hat_minus,bd_F_hat_plas,bd_F_hat_minus
end

function setup_RHS(cell_Enum,cell_Fnum,E_hat,F_hat)
    RHS=zeros(size(cell_Enum)[1],4)
    
    for i in 1:size(cell_Enum)[1]
        for j in 1:4
            RHS[i,j] = -(E_hat[cell_Enum[i,3],j]-E_hat[cell_Enum[i,2],j]+F_hat[cell_Fnum[i,3],j]-F_hat[cell_Fnum[i,2],j])
        end
    end
    return RHS
end

function setup_Qcon_hat(Qcon,volume)
    Qcon_hat=ones(size(Qcon)[1],size(Qcon)[2])
    for i in 1:size(Qcon)[1]
        for j in 1:size(Qcon)[2]
            Qcon_hat[i,j]=Qcon[i,j]*volume[i]
        end
    end
    return Qcon_hat
end

function Qhat_to_Q(Qcon_hat,volume)
    Qcon=ones(size(Qcon_hat)[1],size(Qcon_hat)[2])
    for i in 1:size(Qcon_hat)[1]
        for j in 1:size(Qcon_hat)[2]
            Qcon[i,j]=Qcon_hat[i,j]/volume[i]
        end
    end
    return Qcon
end