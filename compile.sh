mkdir -p bin
find . -name "*.java" > source.txt;
javac -cp bin:lib/*:../FlexSC/lib/*:../FlexSC/bin -d bin @source.txt;
