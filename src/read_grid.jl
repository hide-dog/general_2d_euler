# ----------------------
# -- read             --
# ----------------------

function read_nodes(skipnum)
    """ nodes[i][j]
        i:接点番号
        j=0 : 点の番号
        j=1 : 点のx座標
        j=2 : 点のy座標
        j=3 : 点のz座標
    """
    fff=[]
    open("grid/nodes", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_nodes=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    nodes=ones(num_nodes,length(split(fff[1+skipnum]," ")))
    for i in 1:num_nodes
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            nodes[i,j]=parse(Float64,temp[j]) 
        end
    end
    return nodes
end 

function read_elements(skipnum)
    """ elements[i][j](Int)
        i:cell番号
        j=1:接点番号1
        j=2:接点番号2
        j=3:接点番号3
        j=4:接点番号4
        ループするようになっている
    """
    fff=[]
    open("grid/elements", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_cell=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    elements=ones(Int64,num_cell,length(split(fff[1+skipnum]," ")))
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            elements[i,j]=parse(Int64,temp[j]) 
        end
    end
    return elements,num_cell
end 

function read_elements_boundary(skipnum)
    """ elements_bdx[i][j] (Int)
        i:bd番号
        j=0:bd番号
        j=1:接点番号1
        j=2:接点番号2
        j=3:セル番号1
        j=4:セル番号2
    """
    fff=[]
    open("grid/elements_bdx", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_cell=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    elements_bdx=ones(Int64,num_cell,length(split(fff[1+skipnum]," ")))
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            elements_bdx[i,j]=parse(Int64,temp[j]) 
        end
    end

    fff=[]
    open("grid/elements_bdy", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_cell=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    elements_bdy=ones(Int64,num_cell,length(split(fff[1+skipnum]," ")))
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            elements_bdy[i,j]=parse(Int64,temp[j]) 
        end
    end
    return elements_bdx,elements_bdy
end

function read_vecA(skipnum)
    """ vecAx[i][j]
        i:Ax番号
        j=0 : nodes num1 (面を作る点の番号)
        j=1 : nodes num2 (面を作る点の番号)
        j=2 : x vec
        j=3 : y vec
        j=4 : z vec
    """
    fff=[]
    open("grid/vecAx", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_cell=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    vecAx=ones(num_cell,length(split(fff[1+skipnum]," ")))
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            vecAx[i,j]=parse(Float64,temp[j]) 
        end
    end

    fff=[]
    open("grid/vecAy", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_cell=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    vecAy=ones(num_cell,length(split(fff[1+skipnum]," ")))
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            vecAy[i,j]=parse(Float64,temp[j]) 
        end
    end
    return vecAx, vecAy
end 

function read_cell_EFnum(skipnum)
    """ cell_Enum[i][j] (Int)
        i:セル番号
        j=0:セル番号
        j=1:bd番号1
        j=2:bd番号2
    """
    fff=[]
    open("grid/cell_Enum", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_cell=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    cell_Enum=ones(Int64,num_cell,length(split(fff[1+skipnum]," ")))
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            cell_Enum[i,j]=parse(Int64,temp[j]) 
        end
    end

    fff=[]
    open("grid/cell_Fnum", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_cell=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    cell_Fnum=ones(Int64,num_cell,length(split(fff[1+skipnum]," ")))
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            cell_Fnum[i,j]=parse(Int64,temp[j]) 
        end
    end
    return cell_Enum,cell_Fnum
end

function read_volume(skipnum)
    fff=[]
    open("grid/volume", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_cell=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    volume=ones(num_cell,length(split(fff[1+skipnum]," ")))
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            volume[i,j]=parse(Float64,temp[j]) 
        end
    end
    return volume
    
end

function read_flux_volume(skipnum)
    
    fff=[]
    open("grid/fluxx_volume", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_cell=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    fx_volume=ones(num_cell,length(split(fff[1+skipnum]," ")))
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            fx_volume[i,j]=parse(Float64,temp[j]) 
        end
    end

    fff=[]
    open("grid/fluxy_volume", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_cell=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    fy_volume=ones(num_cell,length(split(fff[1+skipnum]," ")))
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            fy_volume[i,j]=parse(Float64,temp[j]) 
        end
    end
    
    return fx_volume,fy_volume
end

function read_cell_vecA(skipnum)
    fff=[]
    open("grid/cell_vecAx", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_cell=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    cell_vecAx=ones(num_cell,length(split(fff[1+skipnum]," ")))
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            cell_vecAx[i,j]=parse(Float64,temp[j]) 
        end
    end

    fff=[]
    open("grid/cell_vecAy", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_cell=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    cell_vecAy=ones(num_cell,length(split(fff[1+skipnum]," ")))
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")
        for j in 1:length(temp)
            cell_vecAy[i,j]=parse(Float64,temp[j]) 
        end
    end
    return cell_vecAx,cell_vecAy
end

function read_allgrid()
    skip=1
    nodes=read_nodes(skip)
    elements, num_cell=read_elements(skip)
    ele_bdx,ele_bdy=read_elements_boundary(skip)
    vecAx, vecAy=read_vecA(skip)
    cell_Enum,cell_Fnum=read_cell_EFnum(skip)
    volume=read_volume(skip)
    fx_volume,fy_volume=read_flux_volume(skip)
    cell_vecAx,cell_vecAy=read_cell_vecA(skip)
    println("fin read grid")
    return nodes,elements,ele_bdx,ele_bdy,vecAx,vecAy,num_cell,cell_Enum,cell_Fnum,volume,fx_volume,fy_volume,cell_vecAx,cell_vecAy
end
