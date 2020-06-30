# ----------------------
# -- read             --
# ----------------------
function read_nodenum(skipnum)
    """ 
    xmax : 仮想セルも含めたnodeのxの数
    ymax : 仮想セルも含めたnodeのyの数
    """
    fff=[]
    open("grid/nodesnum", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_nodes=length(fff)-skipnum

    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    temp = split(fff[2]," ")
    xmax = parse(Int64,temp[1]) 
    ymax = parse(Int64,temp[2]) 
    
    return xmax, ymax
end

function read_nodes(skipnum,xmax,ymax)
    """ 
    nodes[i][j][k]
    i : x点の番号
    j : y点の番号
    k=1 : 点のx座標
    k=2 : 点のy座標
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

    nodes=zeros(xmax,ymax,2)
    for i in 1:num_nodes
        temp=split(fff[i+skipnum]," ")

        x = parse(Int64,temp[1])
        y = parse(Int64,temp[2])
        nodes[x,y,1]=parse(Float64,temp[3])
        nodes[x,y,2]=parse(Float64,temp[4]) 
    end
    return nodes
end 

function read_vecA(skipnum,xmax,ymax)
    """ vecAx[i][j][k]
        i : x番号
        j : y番号
        k=1 : x vec
        k=2 : y vec
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

    vecAx=zeros(xmax,ymax-1,2)
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")

        x = parse(Int64,temp[1])
        y = parse(Int64,temp[2])
        vecAx[x,y,1]=parse(Float64,temp[3])
        vecAx[x,y,2]=parse(Float64,temp[4]) 
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

    vecAy=zeros(xmax-1,ymax,2)
    for i in 1:num_cell
        temp=split(fff[i+skipnum]," ")

        x = parse(Int64,temp[1])
        y = parse(Int64,temp[2])
        vecAy[x,y,1]=parse(Float64,temp[3])
        vecAy[x,y,2]=parse(Float64,temp[4]) 
    end
    return vecAx, vecAy
end 

function read_allgrid()
    skip=1
    xmax,ymax=read_nodenum(skip)
    nodes=read_nodes(skip,xmax,ymax)
    vecAx, vecAy=read_vecA(skip,xmax,ymax)
    println("fin read grid")
    return xmax,ymax,nodes,vecAx,vecAy
end
