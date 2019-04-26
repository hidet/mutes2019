# check the mount dir
#mountdir="/s"
#if [ -z "$(ls -A ${mountdir})" ]; then
#    echo "mounting ${mountdir}, need password of oper at sutera..."
#    sshfs oper@sutera.intra.j-parc.jp:/raid/e62 ${mountdir} -o ro
#else
#   echo "${mountdir} exists! go ahead..."
#fi

LOCALDIR="/s/TMU_2019G/"
LNDIR="/Users/tatsuno/work/mutes/data/TMU_2019G/"
#mkdir -p $LNDIR

if [ -d "${LNDIR}" ]; then
    # symbolic link
    echo "creating symbolic links"
    #lndir ${LOCALDIR} ${LNDIR}
    lndir -silent ${LOCALDIR} ${LNDIR}
else
    echo "${LNDIR} does not exit. something is wrong..."
    exit 1
fi
