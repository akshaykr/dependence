declare -a categs="apple car cow cup dog pear tomato horse"
declare -a orients="000-000 045-000 045-090 045-180 066-117 066-333 066-297 090-090 090-158 090-225"

rm -rf tomatoes

for categ in $categs #"@{categs[@]}"
do
  rm -rf $categ"s"
  mkdir $categ"s"
done

for categ in $categs #"@{categs[@]}"
do
  for orient in $orients #"@orients[@]}"
  do
    src="../images/"$categ"s/"$categ"*"$orient".png"
    dest=$categ"s/"
    cp $src $dest
#     echo $src echo $dest
  done
done

mv tomatos tomatoes

