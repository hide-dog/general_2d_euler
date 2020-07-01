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