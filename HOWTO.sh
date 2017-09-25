################################################################################
# this directory contains the source code and the documentation that we provide
# as companion material for the publication
################################################################################

mkdir KAPAC
mkdir PAQR

touch KAPAC/README.md
touch PAQR/README.md

git init
git remote add origin ssh://git@git.scicore.unibas.ch:2222/zavolan_public/PAQR_KAPAC.git
git add .
git commit -m "Initial commit"
git push -u origin master

