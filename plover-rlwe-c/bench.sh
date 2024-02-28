#!/bin/bash

#	Disable frequency scaling until the next boot. Intel:
#		echo 1 > /sys/devices/system/cpu/intel_pstate/no_turbo
#	AMD:
#		echo 0 > /sys/devices/system/cpu/cpufreq/boost

for dut in \
	PLOVER_128_1	PLOVER_128_2	PLOVER_128_4	\
	PLOVER_128_8	PLOVER_128_16	PLOVER_128_32
do
	logf=bench_$dut.txt
	echo === $logf ===
	make obj-clean
	make RACCF="-D"$dut" -DBENCH_TIMEOUT=2.0"
	./xtest | tee $logf
	grep -e '_core_keygen' -e '_core_sign' -e '_core_verify' plover_core.su >> $logf
done
