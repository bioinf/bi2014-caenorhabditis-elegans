1. Готовим выборку из белков (для получения более менее быстрого результата желательно не больше 200)
2. Прогнать её через exonerate с параметрами:
exonerate --model protein2genome --query FILE_WITH_PROTEINS --target FILE_WITH_GENOME --bestn 1 --percent 80 --showvulgar no --showalignment no --showtargetgff yes > FILE_WITH_ALIGNMENT.gff
3. Идём в папку augustus/scripts, в ней удаляем файл gff2gbSmallDNA.pl и заменяем его на файлик с таким же именем, который лежит в папке AUGUSTUS на dropbox
4. Перегоняем наш результат от exonerate из gff в genbank:
./gff2gbSmallDNA.pl FILE_WITH_ALIGNMENT.gff FILE_WITH_GENOME 1000 ./genes.gb

*я к сожалению так и не поняла, что значит 3ий параметр, т.е. 1000 =) но задавала как на большинстве сайтов

5. Следуем инструкции с http://www.molecularevolution.org/molevolfiles/exercises/augustus/training.html
(только команде randomSplit.pl genes.gb 100 - нужно число указать так, чтобы в train оказалось >= 100)


