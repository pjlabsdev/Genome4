# Seed-Suche: längste exakte Treffer (SMEM)

Stand: gebaut und geprüft (`gc-index`, `smem.rs`). Messung:
[M3 Seed-Suche](../../bench/M3-seed-suche.md). Noch nicht gebaut: Nachsuche in langen Treffern
(Re-Seeding), Verketten der Seeds und das Alignment selbst.

## Aufgabe
Ein Read wird nicht als Ganzes gesucht, weil er Lesefehler und echte Varianten enthält. Gesucht
werden die Abschnitte des Reads, die exakt in der Referenz vorkommen und sich in keine Richtung
verlängern lassen, ohne dass der Treffer verloren geht. Solche Abschnitte heißen maximale exakte
Treffer (MEM). Ein SMEM (supermaximal) ist ein MEM, der in keinem anderen MEM desselben Reads
enthalten ist. Die SMEMs sind die Seeds, von denen das Alignment später ausgeht (Li 2013).

## Suche in beide Richtungen mit einer BWT
Die gewöhnliche Rückwärtssuche im [FM-Index](fm-index.md) verlängert ein Muster nur nach links.
Für SMEMs muss man auch nach rechts verlängern. Dafür nutzt Genome 4.0 dieselbe Eigenschaft wie
BWA (FMD-Index, Li 2012): Der Indextext besteht aus dem Vorwärtsstrang und seinem reversen
Komplement und ist damit sein eigenes reverses Komplement. Kommt ein Muster `P` an Textposition
`i` vor, dann kommt `rk(P)` an Position `2L - i - |P|` vor.

Geführt wird deshalb ein Doppelbereich `(k, kr, n)`: `k` ist die erste Zeile der Suffixe, die
mit `P` beginnen, `kr` die erste Zeile für `rk(P)`, und `n` die gemeinsame Zeilenzahl.

- **Nach links um Base c** (aus `P` wird `cP`): `k` und `n` wie in der Rückwärtssuche über
  `occ`. Für `kr` gilt: Im Bereich von `rk(P)` stehen die Fortsetzungen in der Reihenfolge
  Endzeichen, A, C, G, T, und `rk(P)·x` kommt genau so oft vor wie `komp(x)·P`. Man berechnet
  also die Größen aller vier Linksverlängerungen und zählt in `kr` jene davor ab, deren
  Komplement kleiner ist. Das Endzeichen folgt höchstens einmal; seine Zahl ist `n` minus die
  Summe der vier Größen.
- **Nach rechts um Base c** (aus `P` wird `Pc`): Linksverlängerung von `rk(P)` um `komp(c)`,
  mit vertauschten Rollen von `k` und `kr`.

## Ablauf je Read
Wie in BWA-MEM (Li 2013, Funktion `bwt_smem1` in bwa):

1. Ab Read-Position `x` nach rechts verlängern. Jedes Mal, bevor der Bereich kleiner wird, den
   bisherigen Bereich mit seinem Ende merken. Schluss, wenn der Bereich leer wird, ein N kommt
   oder der Read endet.
2. Die gemerkten Bereiche, die am weitesten reichenden zuerst, gemeinsam Schritt für Schritt nach
   links verlängern. Fällt ein Bereich weg, ist er ein SMEM, wenn in diesem Schritt kein weiter
   reichender Bereich überlebt hat und noch keiner mit demselben Anfang gemeldet wurde.
3. Weiter bei der Position, an der der am weitesten reichende Treffer endet.

Damit findet man alle SMEMs des Reads, jeden genau einmal. Ein N im Read beendet jeden Treffer.
`smems_at` nimmt zusätzlich eine Mindesthäufigkeit entgegen; damit lässt sich später die
Nachsuche von BWA-MEM (längere Treffer in der Mitte neu ansetzen und nur häufigere Treffer
zulassen) ohne neue Suche bauen.

Gemeldet werden Treffer ab einer Mindestlänge, Voreinstellung 19 wie in BWA-MEM (`-k 19`).

## Tabelle für kurze Muster
Die ersten Schritte jeder Suche laufen auf sehr breiten Bereichen: Die beiden Enden liegen in
weit entfernten Zählblöcken, jeder Schritt kostet zwei volle Zählungen und zwei Zugriffe auf
fernen Speicher. Am ganzen Genom mit echten Reads sind das 42 von 100 Schritten (Bereich breiter
als 1.000 Zeilen, gezählt an 50.000 Paaren).

Für kurze Muster gibt es nur wenige Möglichkeiten, bei Länge n genau 4^n. Beim Laden des Index
werden deshalb die Bereiche aller Muster bis Länge 10 einmal berechnet, mit genau denselben
Schritten, die sonst die Suche ginge (1.398.100 Einträge zu 24 Bytes, rund 34 MB). Die Suche
schlägt nach, solange das Muster höchstens 10 Basen hat, in beiden Richtungen, und rechnet erst
danach. Bei kleinen Referenzen ist die Tabelle kürzer: nie mehr Muster der längsten Länge als
Zeilen im Index. Die Tabelle steht nicht in der Indexdatei, das Format bleibt gleich.

Am Ergebnis ändert das nichts. Ein Test vergleicht jeden Eintrag mit der einfachen Suche und die
Suche mit Tabelle gegen die Suche ohne Tabelle (`short_table_changes_nothing`).

## Speicherzugriffe vorab anstoßen
Beim Schritt nach links werden alle noch lebenden Bereiche um dieselbe Base verlängert, und diese
Schritte hängen nicht voneinander ab. Am ganzen Genom liegt jeder Zählblock fast immer nicht im
Cache, der Prozessor wartet also bei jedem Schritt auf den Hauptspeicher. Deshalb werden zuerst
die Zählblöcke aller Bereiche dieses Schritts angefordert (`prefetch_row`, über `gc-sys`) und erst
dann gezählt. Die Wartezeiten überlappen sich. Am Ergebnis ändert das nichts, es ist nur ein
Hinweis an den Prozessor. Muster, die in der Tabelle für kurze Muster stehen, brauchen keinen
Block und werden übersprungen. Innerhalb eines Reads hängt der Schritt nach rechts vom vorigen
ab.

**Mehrere Reads verschränkt** (wie BWA-MEM2). Die Suchen sind deshalb in Schritte zerlegt
(`SmemSearch` für die ersten beiden Runden, `ForwardSeed` für die dritte): Jeder Schritt rechnet
bis einschließlich des nächsten Zugriffs auf einen Zählblock. `collect_intervals_many` führt die
drei Runden für 32 Reads zugleich: erst fordert jede Suche die Blöcke ihres nächsten Schritts an,
dann rechnet jede ihren Schritt. `smems_at` und `seed_forward` gehen dieselben Schritte am Stück,
die Ergebnisse sind gleich (Test `many_reads_at_once_give_the_same_intervals`, dazu am ganzen
GRCh38 dieselbe BAM).

**Befund.** Ein Profil (callgrind, chr20) zeigt: Die Seed-Suche wartet weniger auf den Speicher
als angenommen, sie rechnet. Am chr20-Index, der zum Teil im Cache liegt, braucht sie für 200.000
Reads 4,1 s, am ganzen GRCh38 6,4 s auf einem Kern; rund die Hälfte der Befehle steckt im Zählen
(`occ4`). Das Verschränken bringt am ganzen Genom deshalb nur wenige Prozent. Seither zählt
`occ4` in der hinteren Hälfte eines Blocks vom Zählstand des nächsten Blocks zurück, statt vom
Anfang vorwärts, so sind es höchstens 64 statt 128 Basen. Zusammen: Seed-Suche am ganzen GRCh38
6,0 s statt 6,4 s (jeweils bester von drei Läufen auf einem Kern, siehe
[M10](../../bench/M10-ueberlappen.md)). Im ganzen `fq2bam` am ganzen GRCh38 auf 4 Kernen der
Cloud war das allein nicht messbar (Ausrichten 19,7 bis 20,0 s vorher wie nachher).

Seit Zweig `m18-schnellweg` läuft das Verschränken zusammen mit dem
[Schnellweg](schnellweg.md): In jeder Gruppe von 32 Reads gehen zuerst die Reads mit
eindeutigem Treffer über den Schnellweg, nur die übrigen kommen in die verschränkte Suche.
Schalter zurück zum alten Verhalten: `SeedOptions::interleave = false` (dann sucht jeder Read
einzeln). Bei `--schnellweg pruefen` sucht jeder Read immer einzeln. Messung Cloud, 421.564
HG002-Paare am ganzen GRCh38, 4 Kerne, abwechselnd: verschränkt 22,1 / 22,0 / 21,4 s, einzeln
23,4 / 22,5 / 22,4 s, also rund 4 % weniger; die BAM ist byte-gleich.

## Was ein Treffer enthält
Anfang und Ende im Read und den Doppelbereich. Die Zeilenzahl ist die rohe Häufigkeit im
Indextext beider Stränge. Sie kann Fundstellen enthalten, die über eine Contig-Grenze, die
Stranggrenze oder eine gefüllte N-Strecke reichen; diese fallen erst beim Orten weg
(`hit_at`, siehe [FM-Index](fm-index.md)).

## Validierung

| Prüfung | Ergebnis |
|---|---|
| Doppelbereiche gegen die einfache Rückwärtssuche: bei jedem Verlängerungsschritt nach links und rechts müssen `k`, `kr` und `n` den Bereichen von Muster und reversem Komplement gleichen (300 Muster, Zufallsgenom mit Wiederholungen, zwei Contigs, N-Strecke) | gleich |
| `smems_at` gegen eine unabhängige Vergleichssuche, die alle Abschnitte über `x` durch Abzählen im Text bewertet, mit Mindesthäufigkeit 1, 2 und 3 (120 Reads mit Fehlern und N) | gleich |
| SMEMs ganzer Reads gegen die Vergleichssuche über die längste Übereinstimmung ab jeder Position (150 Reads) | gleich |
| chr20: 200.000 erzeugte Reads zu 150 Basen mit 1 % Fehlern, jeder der 471.822 Treffer nachgeprüft: Bereich gleich der einfachen Suche, links und rechts nicht verlängerbar, keiner in einem anderen enthalten | 0 Abweichungen, beide Varianten |

## Grenzen
- Noch kein Vergleich der Seeds mit denen von bwa Read für Read. Der Vergleich mit bwa folgt auf
  der Ebene der fertigen Alignments.
- Noch nicht geprüft: ob bei erzeugten Reads der wahre Ort unter den Orten der Seeds ist.
- Je Verlängerung ein Zugriff auf den Index für beide Enden des Bereichs (`occ4_pair`), ohne
  Vorabruf. Die Suchen mehrerer Reads abwechselnd weiterzuführen, damit der Prozessor auf mehrere
  Stellen des Arbeitsspeichers zugleich wartet, wurde gebaut und gemessen und brachte auf dem
  Messrechner nichts: Die Suche hängt dort am Rechnen, nicht am Speicher. Tempo:
  [M3 Tempo](../../bench/M3-tempo.md), [M6 Tempo](../../bench/M6-tempo.md).

## Literatur
Li 2012 (FMD-Index), Li 2013 (BWA-MEM, SMEM als Seeds), siehe [Literatur](../literatur.md).
