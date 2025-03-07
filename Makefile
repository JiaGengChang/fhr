WDIR=/home/users/nus/e1083772/fhr
PBSOPTS=-o $(WDIR)/.pbs -j oe -l mem=4GB -l ncpus=2 -lwalltime=02:00:01
SCOREDIR ?= "./scores/models/sPLSDA/RNA"

one:
	qsub ${PBSOPTS} -J 49-50:2 ${WDIR}/submit.sh

all:
	qsub ${PBSOPTS} -J 0-49 ${WDIR}/submit.sh

score:
	R -q -e "source('./utils/score_cv.R'); scoredir <- '${SCOREDIR}'; result <- score_repeated_cv(scoredir, n_shuffle=10, n_fold=5); summary(result);"

pbs:
	find ${WDIR}/.pbs/*.OU -type f -delete