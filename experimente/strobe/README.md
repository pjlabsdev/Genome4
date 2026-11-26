# Messwerkzeug Strobemer-Seeds (Experiment, nicht Teil von genome)

Prüft, ob Randstrobe-Seeds (Syncmere k=20, s=16, zweites Strobe 2 bis 9 Syncmere weiter, wie
strobealign, Sahlin 2022) Reads schnell und sicher verorten. Je Read: Strobes beider Stränge im
Index nachschlagen (Strobes mit mehr als `MAXOCC` Orten, Vorgabe 1000, werden gesperrt), Treffer
nach Diagonale gruppieren, die 16 größten Gruppen ohne Lücken (Verschiebung ±3) Base für Base
prüfen.

    cargo build --release
    zcat ref.fa.gz | ./target/release/strobe - reads_1.fq [reads_2.fq]

Ausgabe je FASTQ `strobeN.tsv`: Name, Contig, Ort (1-basiert), Strang (1 = Gegenstrang), Fehlzahl,
Fehlzahl des Zweitbesten, gesperrte Strobes, Treffer der größten Gruppe. Steht die Wahrheit im
Namen (`chr:pos:…`), zählt das Programm richtig und falsch mit. Index ganzes GRCh38: 586,6 Mio. Strobes, 4,8 GB, Spitze 10,2 GB beim Bau.

Messung Cloud (4 Kerne, 300.000 simulierte Reads über GRCh38): 34 µs je Read; 82,4 % „sicher“
(≤ 5 Fehler, Zweitbester ≥ 3 schlechter, nichts gesperrt), davon 7 falsch.
