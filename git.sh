git remote add upstream https://github.com/DelinLi/Bioinfo_Delin.git
git remote -v

git fetch upstream
git checkout master
git merge upstream/master

git add .
git commit -m "upadte"
git push origin master
