#!/bin/bash
# init file for matilda
#
# chkconfig: - 98 98
# description: matilda
#
# processname: matilda

SHELL_SCRIPT_NAME=${BASH_SOURCE:-${0}}

PROJECT_DIR=$(dirname $(readlink -f "${SHELL_SCRIPT_NAME}"))
MANAGE="${PROJECT_DIR}/matilda-manage.sh"
LOGFILE="${PROJECT_DIR}/matilda-manage.log"
PIDFILE="${PROJECT_DIR}/matilda-manage.pid"
EXECUTABLE_SCRIPT="${PROJECT_DIR}/in-screen.sh"
# STARTER_SCRIPT is limited to 15 characters as that is what the script in check_pid_is_running() gets from the check. 
STARTER_SCRIPT=st-matilda.sh
RETVAL=0
SLEEP_DELAY=1.5  # wait for process, sometimes
MATILDA_CONDA_ENV=matilda

REQUIRED_HOST=usaxscontrol.xray.aps.anl.gov

if [ "${HOSTNAME}" != "${REQUIRED_HOST}" ]; then
    echo "Do not run this on ${HOSTNAME}, run it on '${REQUIRED_HOST}' as user 'usaxs'"
    exit 1
fi



activate_conda(){
    if [ "${CONDA_EXE}" == "" ]; then
        echo "Need CONDA_EXE defined to activate '${MATILDA_CONDA_ENV}' environment."
        echo "That is defined by activating *any* conda environment."
        exit 1
    fi
    source ${CONDA_PREFIX:-"/APSshare/miniconda/x86_64"}/etc/profile.d/conda.sh
    CONDA_ENV=matilda
    conda activate "${CONDA_ENV}"      
    #CONDA_ROOT=$(dirname $(dirname $(readlink -f "${CONDA_EXE}")))  
    #source ${CONDA_PREFIX:-"/APSshare/miniconda/x86_64"}/etc/profile.d/conda.sh
    #source "${CONDA_ROOT:-"/APSshare/miniconda/x86_64"}/etc/profile.d/conda.sh"
    #conda activate "${MATILDA_CONDA_ENV}"
}


get_pid(){
    PID=$(/bin/cat "${PIDFILE}")
    return $PID
}


function pid_is_running(){
	get_pid
	if [ "${PID}" == "" ]; then
		# no PID in the PIDFILE
		RETVAL=1
	else
		RESPONSE=$(ps -p ${PID} -o comm=)
		if [ "${RESPONSE}" == "${STARTER_SCRIPT}" ]; then
			# PID matches the Matilda server profile
			RETVAL=0
		else
			# PID is not Matilda server
			RETVAL=1
		fi
	fi
	return "${RETVAL}"
}


start(){
    #activate_conda # done by st-matilda.sh script anyway
    cd "${PROJECT_DIR}"
    "${EXECUTABLE_SCRIPT}" 2>&1 >> "${LOGFILE}" &
    sleep "${SLEEP_DELAY}"
    PID=$(pidof -x ${STARTER_SCRIPT})
    /bin/echo "${PID}" > "${PIDFILE}"
    /bin/echo \
        "# [$(/bin/date -Is) $0] started ${PID}: ${EXECUTABLE_SCRIPT}" \
        2>&1 \
        >> "${LOGFILE}" &
    sleep "${SLEEP_DELAY}"
    tail -1 "${LOGFILE}"
}


stop(){
    get_pid

    if pid_is_running; then
    	/bin/echo "# [$(/bin/date -Is) $0] stopping ${PID}: ${EXECUTABLE_SCRIPT}" 2>&1 >> ${LOGFILE} &
    	kill "${PID}"
    else
		/bin/echo "# [$(/bin/date -Is) $0] not running ${PID}: ${EXECUTABLE_SCRIPT}" 2>&1 >> ${LOGFILE} &
    fi
    sleep "${SLEEP_DELAY}"
    tail -1 "${LOGFILE}"

    /bin/cp -f /dev/null "${PIDFILE}"
}


restart(){
    stop
    start
}


status(){
    if pid_is_running; then
		echo "# [$(/bin/date -Is) $0] running fine, so it seems"
    else
		echo "# [$(/bin/date -Is) $0] could not identify running process ${PID}"
    fi
}


checkup(){
    # 'crontab -e` to add entries for automated (re)start
    #=====================
    # call periodically (every 5 minutes) to see if Matilda server is running
    #=====================
    #	     field	    allowed values
    #	   -----	  --------------
    #	   minute	  0-59
    #	   hour 	  0-23
    #	   day of month   1-31
    #	   month	  1-12 (or names, see below)
    #	   day of week    0-7 (0 or 7 is Sun, or use names)
    #
    # */5 * * * * /home/beams/JEMIAN/Documents/projects/BCDA-APS/tiled-template/tiled-manage.sh checkup 2>&1 > /dev/null

    if pid_is_running; then
		echo "# [$(/bin/date -Is) $0] running fine, so it seems" 2>&1 > /dev/null
    else
		echo "# [$(/bin/date -Is) $0] could not identify running process ${PID}, starting new process" 2>&1 >> "${LOGFILE}"
		start
    fi
}


case "$1" in
    start)    start ;;
    stop)     stop ;;
    restart)  restart ;;
    checkup)  checkup ;;
    status)   status ;;
    *)
        echo $"Usage: $0 {start|stop|restart|checkup|status}"
        exit 1
esac
