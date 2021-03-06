function set_boundary(Qbase,cellxmax,cellymax,vecAx,vecAy,bdcon)
    """bdcon[i][j]
        i : 境界番号(x-,x+ ,y-,y+)
        j=1-6 : "bd1_con":"2",
                "bd1_rho":"1.0",
                "bd1_u":"300.0",
                "bd1_v":"0.0",
                "bd1_p":"1.0",
                "bd1_T":"300.0",
    """
    # bd1 = x-
    if Int(bdcon[1][1]) == 0
        for j in 1:cellymax
            for k in 1:4
                Qbase[1,j,k] = bdcon[1][k+1]
            end
        end
    elseif Int(bdcon[1][1]) == 1
        for j in 1:cellymax
            for k in 1:4
                Qbase[1,j,k] = Qbase[2,j,k]
            end
        end
    elseif Int(bdcon[1][1]) == 2
        for j in 1:cellymax
            xvec = vecAx[2,j,1]
            yvec = vecAx[2,j,2]
            u = Qbase[2,j,2]
            v = Qbase[2,j,3]

            Qbase[1,j,1] = Qbase[2,j,1]
            Qbase[1,j,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[1,j,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
            Qbase[1,j,4] = Qbase[2,j,4]  
        end
    elseif Int(bdcon[1][1]) == 3
        for j in 1:cellymax
            u = Qbase[2,j,2]
            v = Qbase[2,j,3]

            Qbase[1,j,1] = Qbase[2,j,1]
            Qbase[1,j,2] = -u
            Qbase[1,j,3] = -v
            Qbase[1,j,4] = Qbase[2,j,4]  
        end
    else
        println("\n check boundary condition ! \n")
        throw(UndefVarError(:x))
    end

    # bd2 = x+
    if Int(bdcon[2][1]) == 0
        for j in 1:cellymax
            for k in 1:4
                Qbase[cellxmax,j,k] = bdcon[2][k+1]
            end
        end
    elseif Int(bdcon[2][1]) == 1
        for j in 1:cellymax
            for k in 1:4
                Qbase[cellxmax,j,k] = Qbase[cellxmax-1,j,k]
            end
        end
    elseif Int(bdcon[2][1]) == 2
        for j in 1:cellymax
            xvec = vecAx[cellymax,j,1]
            yvec = vecAx[cellymax,j,2]
            u = Qbase[cellxmax-1,j,2]
            v = Qbase[cellxmax-1,j,3]

            Qbase[cellxmax,j,1] = Qbase[cellxmax-1,j,1]
            Qbase[cellxmax,j,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[cellxmax,j,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
            Qbase[cellxmax,j,4] = Qbase[cellxmax-1,j,4]  
        end
    elseif Int(bdcon[2][1]) == 3
        for j in 1:cellymax
            u = Qbase[cellxmax-1,j,2]
            v = Qbase[cellxmax-1,j,3]

            Qbase[cellxmax,j,1] = Qbase[cellxmax-1,j,1]
            Qbase[cellxmax,j,2] = -u
            Qbase[cellxmax,j,3] = -v
            Qbase[cellxmax,j,4] = Qbase[cellxmax-1,j,4]  
        end
    else
        println("\n check boundary condition ! \n")
        throw(UndefVarError(:x))
    end

    # bd3 = y-
    if Int(bdcon[3][1]) == 0
        for i in 1:cellxmax
            for k in 1:4
                Qbase[i,1,k] = bdcon[3][k+1]
            end
        end
    elseif Int(bdcon[3][1]) == 1
        for i in 1:cellxmax
            for k in 1:4
                Qbase[i,1,k] = Qbase[i,2,k]
            end
        end
    elseif Int(bdcon[3][1]) == 2
        for i in 1:cellxmax
            xvec = vecAy[i,2,1]
            yvec = vecAy[i,2,2]
            u = Qbase[i,2,2]
            v = Qbase[i,2,3]

            Qbase[i,1,1] = Qbase[i,2,1]
            Qbase[i,1,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[i,1,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
            Qbase[i,1,4] = Qbase[i,2,4]  
        end
    elseif Int(bdcon[3][1]) == 3
        for i in 1:cellxmax
            u = Qbase[i,2,2]
            v = Qbase[i,2,3]

            Qbase[i,1,1] = Qbase[i,2,1]
            Qbase[i,1,2] = -u
            Qbase[i,1,3] = -vecAy
            Qbase[i,1,4] = Qbase[i,2,4]  
        end
    else
        println("\n check boundary condition ! \n")
        throw(UndefVarError(:x))
    end

    # bd4 = y+
    if Int(bdcon[4][1]) == 0
        for i in 1:cellxmax
            for k in 1:4
                Qbase[i,cellymax,k] = bdcon[4][k+1]
            end
        end
    elseif Int(bdcon[4][1]) == 1
        for i in 1:cellxmax
            for k in 1:4
                Qbase[i,cellymax,k] = Qbase[i,cellxmax-1,k]
            end
        end
    elseif Int(bdcon[4][1]) == 2
        for i in 1:cellxmax
            xvec = vecAy[i,cellymax,1]
            yvec = vecAy[i,cellymax,2]
            u = Qbase[i,cellymax-1,2]
            v = Qbase[i,cellymax-1,3]

            Qbase[i,cellymax,1] = Qbase[i,cellymax-1,1]
            Qbase[i,cellymax,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[i,cellymax,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
            Qbase[i,cellymax,4] = Qbase[i,cellymax-1,4]  
        end
    elseif Int(bdcon[4][1]) == 3
        for i in 1:cellxmax
            u = Qbase[i,cellymax-1,2]
            v = Qbase[i,cellymax-1,3]

            Qbase[i,cellymax,1] = Qbase[i,cellymax-1,1]
            Qbase[i,cellymax,2] = -u
            Qbase[i,cellymax,3] = -v
            Qbase[i,cellymax,4] = Qbase[i,cellymax-1,4]  
        end
    else
        println("\n check boundary condition ! \n")
        throw(UndefVarError(:x))
    end

    return Qbase
end