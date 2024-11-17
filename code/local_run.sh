L=100
kappa=0.0
A=1.0
invK=0.0
#for kappa in $(seq 0.0 1.0 20.0)
#do
#    echo "Running L=$L kappa=$kappa, f=$f, g=$g"
#    nohup ./semiflexible_polymer $L $kappa $f $g lo &
#done
kappa=10.0
for L in 10 50 100 200
do
    for kappa in 5.0
    do
        for invK in 1.0 2.0 5.0 10.0 20.0
        do
            echo "Running L=$L kappa=$kappa, A=$A, invK=$invK"
            nohup ./charged_polymer $L $kappa $A $invK 1 ~/Work/charged_polymer/data_hpc/data_pool &
        done
    done
done


L=100
kappa=0.0
A=1.0
invK=0.0
#for kappa in $(seq 0.0 1.0 20.0)
#do
#    echo "Running L=$L kappa=$kappa, f=$f, g=$g"
#    nohup ./semiflexible_polymer $L $kappa $f $g lo &
#done
kappa=10.0
for kappa in 2.0 5.0 10.0
do
    for A in 1.0 5.0
    do
        for invK in 1.0 2.0 5.0 10.0 20.0
        do
            echo "Running L=$L kappa=$kappa, A=$A, invK=$invK"
            nohup ./charged_polymer $L $kappa $A $invK 1 ~/Work/charged_polymer/data_hpc/data_pool &
        done
    done
done