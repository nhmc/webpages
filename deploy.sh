rm -rf deploy
hyde gen
cd ~/Code/Repo/nhmc.github.io
rm -rf *
cp -r ~/Code/Repo/webpages/deploy/* .
git add .
git commit -m'update'
git push
