#! /bin/bash -x

function test_command(){
    $@
    result=$?
    echo $result
    if [ $result -ne 0 ];
    then
        echo "Command $@ failed with exit code $result"
        exit 1
    fi
}



test_command rsat -h
test_command rsat oligo-analysis -h ## Test a subcommand running a Perl script
test_command rsat random-seq -l 100 -n 2 ## Test a perl script
test_command rsat random-seq -l 100 -n 2 | rsat purge-sequence 
test_command rsat random-motif -l 10 ## Test a subcommand running a python script
test_command rsat feature-map -h ## check the GD dependency for feature-map
test_command rsat info-gibbs -h
test_command rsat count-words -h
test_command rsat matrix-scan-quick -h
test_command rsat matrix-clustering -h ## Test if the specific Perl dependencies have been successfully included for this tool
test_command make -f $(pwd)/makefiles/subcommand_tests.mk randseq ## Quick test for one Perl script
test_command make -f $(pwd)/makefiles/subcommand_tests.mk xygraph ## test GD dependency
