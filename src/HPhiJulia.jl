module HPhiJulia
    export HPhi
    currentd = pwd()
    const models = ["Fermion Hubbard", "Spin", "Kondo Lattice", "Fermion HubbardGC","SpinGC", "Kondo LatticeGC", "SpinGCCMA"]
    const methods = ["Lanczos", "TPQ", "Full Diag", "CG", "Time-Evolution"]
    const lattices = ["Chain Lattice", "Square Lattice", "Triangular Lattice", "Honeycomb Lattice", "Ladder", "Kagome"]

       
    cd(homedir()*"/.julia/packages")
    if isdir("HPhiJulia")
    else
        mkdir("HPhiJulia")
    end
    cd("HPhiJulia")

    const HPhipath = homedir()*"/.julia/packages/HPhiJulia/HPhi.build/src/HPhi" 
    const HPhiversion = "3.2.0"
#    const HPhiversion = "3.1.2"
    const HPhitar = "HPhi-"*HPhiversion*".tar.gz"
    const HPhiurl = "https://github.com/issp-center-dev/HPhi/releases/download/v"*HPhiversion*"/"
              
    if isfile(HPhipath) == false
        using CMake
        if isfile(HPhitar) == false
            println("Downloading HPhi...")
            download( HPhiurl*HPhitar,HPhitar)
	        println("done.")
	  
	    end
        run(`tar xzvf $HPhitar`)
        if isdir("HPhi.build") == false
       	    mkdir( "HPhi.build" )
        end
        cd("HPhi.build")
        run(`$cmake ../HPhi-$HPhiversion/`)
        run(`make`)
    end
    cd(currentd)

    function create_stanin(model,lattice,W,L,twoSz,method,
        mu,t,tp,U,V,Vp,J,Jp,nelec)
       fp = open("stan.in","w")
       println(fp,"model = \"$model\"")
       println(fp,"method = \"$method\"")
       println(fp,"lattice = \"$lattice\"")
       if lattice != "Chain Lattice"
              println(fp,"W= $W")
       end
       println(fp,"L = $L")
                     
       if mu != 0 println(fp,"mu = $mu") end
       if t != 0 println(fp,"t = $t") end
       if tp != 0 println(fp,"t` = $tp") end
       if U != 0 println(fp,"U = $U") end
       if V != 0 println(fp,"V = $V") end
       if Vp != 0 println(fp,"V` = $Vp") end
       if J != 0 println(fp,"J = $J") end
       if Jp != 0 println(fp,"J` = $Jp") end
       if nelec != 0 println(fp,"nelec = $nelec")
              if model == "Fermion HubbardGC" || model == "Spin" || model == "SpinGC"
                     error("nelec should not be used in model: $model .")                                   
              end
       end
       if model == "Fermion HubbardGC" || model == "SpinGC"
       else
              println(fp,"2Sz = $twoSz")
       end
       close(fp)
    end

    function read_calcmod()
       calcmod = readlines("calcmod.def")
       Calcmod = Dict(
       "CalcType" => parse(Int8,split(calcmod[5])[2]),
       "CalcModel" => parse(Int8,split(calcmod[6])[2]),
       "ReStart" => parse(Int8,split(calcmod[7])[2]),
       "CalcSpec" => parse(Int8,split(calcmod[8])[2]),
       "CalcEigenVec" => parse(Int8,split(calcmod[9])[2]),
       "InitialVecType" => parse(Int8,split(calcmod[10])[2]),
       "InputEigenVec" => parse(Int8,split(calcmod[11])[2]),
       "OutputEigenVec" => parse(Int8,split(calcmod[12])[2])
       )
       return Calcmod       
    end

    function write_calcmod(Calcmod)
       fp = open("calcmod.def","w")
       println(fp,"#CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG, 4:Time-evolution")
       println(fp,"#CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC")
       println(fp,"#Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart")
       println(fp,"#CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, 4:Restart, 5:Restart&Save")
       println(fp,"CalcType   $(Calcmod["CalcType"])")
       println(fp,"CalcModel   $(Calcmod["CalcModel"])")
       println(fp,"CalcSpec   $(Calcmod["CalcSpec"])")
       println(fp,"CalcEigenVec   $(Calcmod["CalcEigenVec"])")
       println(fp,"InitialVecType   $(Calcmod["InitialVecType"])")
       println(fp,"InputEigenVec   $(Calcmod["InputEigenVec"])")
       println(fp,"OutputEigenVec   $(Calcmod["OutputEigenVec"])")
       close(fp)

       return    
    end


    function HPhi(model,lattice,W,L;expart=false,twoSz=0,method="CG",mpinum = 1,
        mu = 0,t=0,tp=0,U=0,V=0,Vp=0,J=0,Jp=0,nelec=0,OutputEigenVec=false,usenamelist=false)   

        if expart == false
              if  OutputEigenVec
                     error("export should be true when OutputEigenVec = true")
              end
              create_stanin(model,lattice,W,L,twoSz,method,mu,t,tp,U,V,Vp,J,Jp,nelec)             
              if mpinum == 1
                     run(`$HPhipath -s stan.in`)
              else
                     run(`mpirun -np $mpinum $HPhipath -s stan.in`)
              end       
        else
              if usenamelist == false
                     create_stanin(model,lattice,W,L,twoSz,method,mu,t,tp,U,V,Vp,J,Jp,nelec)  
                     run(`$HPhipath -sdry stan.in`)                     
              end
              Calcmod = read_calcmod()
              if  OutputEigenVec                     
                     Calcmod["OutputEigenVec"] = 1
              end 
              write_calcmod(Calcmod)
              if mpinum == 1
                     run(`$HPhipath -e namelist.def`)
              else
                     run(`mpirun -np $mpinum $HPhipath -e namelist.def`)                     
              end
        end

        zvo_energy = readlines("./output/zvo_energy.dat")   
                                       
        Energy = parse(Float64,split(zvo_energy[2])[2])
        Doublon = parse(Float64,split(zvo_energy[3])[2])
        Sz = parse(Float64,split(zvo_energy[4])[2])
        println("Energy is $Energy")
        println("Doublon is $Doublon")
        println("Sz is $Sz")
        results = Dict("Energy" => Energy,"Doublon" => Doublon,"Sz" => Sz)

        if OutputEigenVec
            fname = "output/zvo_eigenvec_0_rank_0.dat" 
            println("loading $fname")
            evec,Hilbert = load_evec(fname)
            return results,evec
       else
              return results
        end

        
        
    end

    function load_evec(fname)
       io = open(fname)
       Lanczos = read(io, Int32)
       Hilbert= read(io, Int64)
       discard= read(io, ComplexF64)
       evec=zeros(ComplexF64, Hilbert)
       for i in 1:Hilbert
              v = read(io, ComplexF64)
              evec[i]= v
       end
       println("Eigenvector was loaded")
       return evec,Hilbert
    end

end