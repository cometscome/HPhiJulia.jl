module HPhiJulia
       #export HPhi
       currentd = pwd()
       
       cd(homedir()*"/.julia/packages")
       if isdir("HPhiJulia")
       else
              mkdir("HPhiJulia")
       end
       cd("HPhiJulia")

       const HPhipath = homedir()*"/.julia/packages/HPhiJulia/HPhi.build/src/HPhi" 
       const HPhiversion = "3.1.2"
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

       function HPhi(model,lattice,W,L;expart=false,twoSz=0,method="CG",mpinum = 1,
       mu = 0,t=0,tp=0,U=0,V=0,Vp=0,J=0,Jp=0,nelec=0)       
              if expart == false
                     
                     fp = open("stan.in","w")
                     println(fp,"model = \"$model\"")
                     println(fp,"method = \"$method\"")
                     println(fp,"lattice = \"$lattice\"")
                     println(fp,"W= $W")
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

              else

              end
              if mpinum == 1
                     run(`$HPhipath -s stan.in`)
              else
                     run(`mpirun -np $mpinum $HPhipath -s stan.in`)
              end
              zvo_energy = readlines("./output/zvo_energy.dat")   
                                       
              Energy = parse(Float64,split(zvo_energy[2])[2])
              Doublon = parse(Float64,split(zvo_energy[3])[2])
              Sz = parse(Float64,split(zvo_energy[4])[2])
              println("Energy is $Energy")
              println("Doublon is $Doublon")
              println("Sz is $Sz")
              results = Dict("Energy" => Energy,"Doublon" => Doublon,"Sz" => Sz)
              return results
       end

end