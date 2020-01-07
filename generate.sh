#!/usr/bin/bash
num_det_blocks=10

#for i in `seq 0 $num_det_blocks`
#do
#  echo "COMPONENT fw_azim${i} = Arm()"
#  echo "AT(0,0,0) RELATIVE mm"
#  echo "ROTATED (0,0,ETA[${i}]) RELATIVE mm"
#  echo "COMPONENT fw_elevation${i} = Arm()"
#  echo "AT(0,0,0) RELATIVE mm"
#  echo "ROTATED (-XI,0,0) RELATIVE fw_azim${i}"
#done
#
#for i in `seq 0 $num_det_blocks`; do
#  echo "COMPONENT bw_azim${i} = Arm()"
#  echo "AT(0,0,0) RELATIVE mm"
#  echo "ROTATED (0,180,ETA[${i}]) RELATIVE mm"
#  echo "COMPONENT bw_elevation${i} = Arm()"
#  echo "AT(0,0,0) RELATIVE mm"
#  echo "ROTATED (-XI,0,0) RELATIVE bw_azim${i}"
#done

for i in `seq 0 $num_det_blocks`; do
 if [ "$i" == "0" ]; then
echo "//COMPONENT fw_tpsd0 = PSD_monitor_TOF(restore_neutron=1,filename=\"fw_tpsd0\",nx=8,ny=2,tmin=0, tmax=10, nt=200, xwidth=0.16, yheight=0.04)"
echo "COMPONENT fw_tpsd0 = TOF_monitor(restore_neutron=1,filename=\"fw_tpsd0\", tmin=0, tmax=10, nt=200,xwidth=0.16, yheight=0.04)"
echo "//COMPONENT fw_tpsd0 = PSD_monitor(restore_neutron=1,filename=\"fw_tpsd0\", xwidth=0.16, yheight=0.04)"
 else
    echo "COMPONENT fw_tpsd${i} = COPY(fw_tpsd0)(filename=\"fw_tpsd${i}\")"
 fi
 echo "AT (0,0,SDD) RELATIVE fw_elevation${i}"
 echo "GROUP d2"
 echo "COMPONENT bw_tpsd${i} = COPY(fw_tpsd0)(filename=\"bw_tpsd${i}\")"
 echo "AT (0,0,SDD) RELATIVE bw_elevation${i}"
 echo "GROUP d2"
done 

