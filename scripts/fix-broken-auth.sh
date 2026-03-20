# files=$(grep -l "cannot generate kerberos ticket in a non-interactive session" working/*/*bilby*/log_data_generation/*.err)

# for file in $files; do
#     job=$(dirname $(dirname $file))
#     event=$(basename $(dirname $job))
#     job=$(basename $job)

#     echo $event $job
    
#     id=$(grep -o -e "[0-9]*.000.000" working/$event/$job/submit/dag_*.submit.dagman.log | tail -n1 | cut -d . -f 1)

#     echo $id
    
#     condor_rm $id
#     sleep 60s
#     rm -rf working/$event/$job
#     asimov production set $event $job -s ready
# done

files=$(ls checkouts/*/C01_offline/*bilby*.ini)
for file in $files; do
    job=$(basename $file | cut -d . -f 1)
    event=$(basename $(dirname $(dirname $file)))
    echo $event $job
    if [ ! -d working/$event/$job ]; then
	asimov production set $event $job -s ready
    fi
    if [ ! -d working/$event/$job/final_result ]; then
	rm -r working/$event/$job/
	asimov production set $event $job -s ready
    fi
done;
