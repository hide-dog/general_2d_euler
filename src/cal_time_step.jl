function time_integration_explicit(dt,Qval_hat,RHS)
    new_Qval_hat=zeros(size(Qval_hat)[1],size(Qval_hat)[2])
    for i in 1:size(Qval_hat)[1]
        for j in 1:size(Qval_hat)[2]
            new_Qval_hat[i,j]=Qval_hat[i,j]+dt*RHS[i,j]
        end
    end
    return new_Qval_hat
end 