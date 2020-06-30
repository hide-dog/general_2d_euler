using Printf

function main()
    output_dir_name="post_result"
    input_result_dir="result"

    #vtk_頭の文字
    front="# vtk DataFile Version 2.0"*"\n"
    back="\n"*"ASCII"*"\n"*"DATASET UNSTRUCTURED_GRID"*"\n"
    write_file(output_dir_name,input_result_dir,front,back)
end
    
function write_file(outdir,inresult_dir,front,back)

    nodes_xmax, nodes_ymax = read_nodenum(1)
    nodes = read_nodes(1)

    elements = read_elements(1)

    inf=readdir(inresult_dir)
    T = zeros(nodes_xmax-3,nodes_ymax-3)
    
    for i in 1:length(inf)
        if occursin(".dat", inf[i]) == true
            fff=[]
            open(inresult_dir*"/"*inf[i], "r") do f
                fff=read(f,String)
            end 
            fff=split(fff,"\n",keepempty=false)

            for i in 2:length(fff)
                fff[i]=replace(fff[i]," \r" => "")
                #temp = split(fff[2]," ")
                T[i-1] = parse(Float64,fff[i])
            end
            dname = replace(inf[i],".dat" => "")
            out_file = replace(inf[i],".dat" => ".vtk")
            
            write_points(nodes,out_file,outdir,dname,front,back)
            write_cells(elements,out_file,outdir)
            write_result(T,out_file,outdir)
        end
    end
end 

function read_nodenum(skipnum)
    # xmax : 仮想セルも含めたnodeのxの数
    # ymax : 仮想セルも含めたnodeのyの数    
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

function read_nodes(skipnum)
    
    fff=[]
    open("grid/nodes_forvtk", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_nodes=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    nodes=zeros(num_nodes,3)
    for i in 1:num_nodes
        temp=split(fff[i+skipnum]," ")

        # x = parse(Float64,temp[1])
        # y = parse(Float64,temp[2])
        # z = parse(Float64,temp[3])
        x = parse(Float64,temp[2])
        y = parse(Float64,temp[3])
        z = 0.0
        nodes[i,1] = x
        nodes[i,2] = y
        nodes[i,3] = z
    end
    return nodes
end 

function read_elements(skipnum)
    fff=[]
    open("grid/element_forvtk", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_elements=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    elements = zeros(num_elements,4)
    for i in 1:num_elements
        temp=split(fff[i+skipnum]," ")

        elements[i,1] = parse(Int64,temp[2])-1
        elements[i,2] = parse(Int64,temp[3])-1
        elements[i,3] = parse(Int64,temp[4])-1
        elements[i,4] = parse(Int64,temp[5])-1
    end
    return elements
end 

function  write_points(nodes,out_file,outdir,dname,front,back)
    a = size(nodes)[1]
    a_st = @sprintf("%1.0f", a)

    fff = outdir*"/"*out_file
    open(fff,"w") do f
        write(f,front)
        write(f,dname)
        write(f,back)
        write(f,"POINTS "*a_st*" float\n")
        for i in 1:a
            x = @sprintf("%7.7f", nodes[i,1])
            y = @sprintf("%7.7f", nodes[i,2])
            z = @sprintf("%7.7f", nodes[i,3])
            write(f, x*" "*y*" "*z*"\n")
        end
    end
    println("fin writing points")
end

function write_cells(elements,out_file,outdir)
    a = size(elements)[1]
    b = a*5
    a_st = @sprintf("%1.0f", a)
    b_st = @sprintf("%1.0f", b)

    fff = outdir*"/"*out_file
    open(fff,"a") do f
        write(f, "CELLS "*a_st*" "*b_st*"\n")
        for i in 1:a
            d1 = @sprintf("%1.0f", elements[i,1])
            d2 = @sprintf("%1.0f", elements[i,2])
            d3 = @sprintf("%1.0f", elements[i,3])
            d4 = @sprintf("%1.0f", elements[i,4])
            write(f, "4 "*d1*" "*d2*" "*d3*" "*d4*"\n")
        end
        write(f, "CELL_TYPES "*a_st*"\n")
        for i in 1:a
            write(f, "9\n")          #四角のみ
        end
    end     
    println("fin writing cells")
end

function write_result(T,out_file,outdir)
    a = length(T)
    a_st = @sprintf("%1.0f", a)

    dtype = "cell_scalars"
    
    temp1 = "CELL_DATA "*a_st*"\n"
    temp2 = "SCALARS "*dtype*" float\nLOOKUP_TABLE default\n"

    fff = outdir*"/"*out_file
    open(fff,"a") do f
        write(f, temp1)
        write(f, temp2)
        for i in 1:a
            data = @sprintf("%7.7f", T[i])
            write(f, data)
            write(f, "\n")
        end
    end 
    print("fin writing "*out_file)
end
# main
main()
