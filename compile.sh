CP="bin:lib/*:../FlexSC/bin:../FlexSC/lib*"

mkdir -p bin
find . -name "*.java"  > source.txt;
javac -cp "$CP" -d bin @source.txt;
rm source.txt
