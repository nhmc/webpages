rm -rf deploy
hyde gen
touch deploy/.nojekyll
cd ~/Code/repo/nhmc.github.io
rm -rf *
cp -r ~/webpages/deploy/* .
git add --all
git commit -m'update'
git push
cd ~/webpages
