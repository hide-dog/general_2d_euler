function set_res(Qcon_hat,Lx,Ly,Ux,Uy,RHS,cellxmax,cellymax)
    res = zeros(cellxmax, cellymax,4)
    
    for i in 2:cellxmax
        for j in 2:cellymax
            Ax_star = zeros(4)
                for l in 1:4
                    for m in 1:4
                        Ax_star[l] = Ax_star[l] + Lx[i-1,j,l,m]*Qcon_hat[i-1,j,m]
                                    + Ly[i,j-1,l,m]*Qcon_hat[i,j-1,m]
                                    + Ux[i+1,j,l,m]*Qcon_hat[i+1,j,m]
                                    + Uy[i,j+1,l,m]*Qcon_hat[i,j+1,m]
                    end
                    res[i,j,l] = Ax_star[l] - RHS[i,j,l]
                end
        end
    end 
    return res
end

function check_converge(res,RHS,cellxmax,cellymax)
    norm2 = zeros(4)

    tempAxb = zeros(4)
    tempb = zeros(4)
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for l in 1:4
                tempAxb[l] = tempAxb[l] + res[i,j,l]^2
                tempb[l] = tempb[l] + RHS[i,j,l]^2
            end
        end
    end

    for l in 1:4
        norm2[l] = (tempAxb[l]/tempb[l])^0.5
    end
    return norm2
end