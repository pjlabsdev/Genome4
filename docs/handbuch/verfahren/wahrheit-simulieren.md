# Wahrheit erzeugen und Ausrichtung prüfen

Stand: gebaut. Werkzeug: [`genome-validate`](../werkzeuge/genome-validate.md).

## Frage
Landet jeder Read dort, wo er herkommt, und stimmt die MAPQ mit der tatsächlichen Fehlerrate
überein? Ein Vergleich mit bwa allein beantwortet das nicht: Gleich wie bwa kann auch gleich falsch
sein. Gebraucht wird eine Wahrheit, die von keinem der beiden Ausrichter abhängt.

## Probe mit zwei Haplotypen
Auf jeden Contig der Referenz werden Varianten gesetzt, Abstände geometrisch mit der Rate
`snv-rate + indel-rate` je Base. Eine Variante ist mit 60 % mischerbig (je zur Hälfte auf
Haplotyp 1 oder 2) und mit 40 % reinerbig. SNVs nehmen eine der drei anderen Basen. Indels haben
die Länge 1 bis `max-indel`, jede weitere Base halb so wahrscheinlich wie die vorige, Lücke oder
Einfügung je zur Hälfte. Orte mit N bekommen keine Variante, und hinter einer Lücke beginnt die
nächste Variante frühestens eine Base nach ihrem Ende. Jeder Contig hat seine eigene Saat, damit
die Probe nicht von der Zahl der Reads abhängt.

## Stück und wahrer Ort
Ein Paar entsteht so: Ort auf dem ganzen Genom gleichverteilt (also Contigs nach Länge
gewichtet), Haplotyp 1 oder 2, Stücklänge normalverteilt (mindestens eine Read-Länge). Ab dem Ort
wird das Stück des Haplotyps gebaut, Base für Base mit ihrem Ort auf der Referenz:

- SNV: andere Base, gleicher Ort.
- Lücke: die Referenzbasen fallen aus, der Ort springt.
- Einfügung: neue Basen vor der Referenzbase, sie tragen deren Ort.
- Beginnt das Stück mitten in einer Lücke seines Haplotyps, beginnt es hinter ihr.

Stücke mit N oder über das Contig-Ende werden verworfen und gezählt. Read 1 und Read 2 sind der
Anfang und das umgedrehte Gegenstück des Endes; welcher von beiden vorne liest, ist zufällig. Das
wahre linke Ende eines Reads ist der Referenzort seiner ersten Base auf dem Vorwärtsstrang.

## Lesefehler
Je Base wird zuerst die Qualität gezogen: höchste Stufe, sonst eine von drei niedrigen (15 %,
30 %, 55 %). Der Anteil der niedrigen steigt quadratisch von `low_start` am Read-Anfang auf
`low_end` am Ende. Danach ist die Base mit genau der Wahrscheinlichkeit 10^(-Q/10) falsch. So
hält jede Qualität, was sie verspricht, wie bei einem gut kalibrierten Gerät. Selten (Illumina
5 auf 100.000 Basen, MGI 2 auf 100.000) fehlt eine Base oder kommt eine hinzu.

| Profil | Stufen | niedrig am Anfang | niedrig am Ende |
|---|---|---|---|
| `illumina` | 2, 12, 23, 37 (wie NovaSeq) | 3 % | 15 % |
| `mgi` | 5, 15, 28, 40 | 2 % | 12 % |

Die Zufallsfolge ist ein fest eingebauter xoshiro256** mit SplitMix64-Start, damit sich die
Daten nie durch ein Update einer fremden Crate ändern.

## Prüfung
Gezählt wird jeder primäre Datensatz (sekundäre und supplementäre nicht). Richtig ist er, wenn
der Contig stimmt und das linke Ende samt abgeschnittener Basen höchstens `toleranz` Basen (20)
neben dem wahren liegt. Abgeschnittene Basen werden zurückgerechnet, damit ein Read mit
Lesefehlern am Anfang nicht als falsch zählt. Die Toleranz deckt Indels der Probe in den ersten
Basen eines Reads.

Je MAPQ stehen die Reads mit Ort, die falschen und die erwarteten falschen: Eine MAPQ von q
verspricht eine Fehlerrate von 10^(-q/10). Ein Ausrichter ist gut kalibriert, wenn falsche und
erwartete in jeder Stufe nahe beieinander liegen. Die Hauptkennzahl aus dem Plan (Abschnitt 0.4)
ist der Anteil falscher Reads unter denen mit MAPQ ab 20.

## Grenzen
Siehe die Werkzeugseite. Das Wichtigste: Eine Probe ohne Strukturvarianten und ohne Sequenz
außerhalb der Referenz ist leichter als ein echter Mensch, und das Fehlermodell ist nicht aus
echten Läufen gelernt. Die Zahlen vergleichen deshalb Ausrichter untereinander auf denselben
Reads; als absolute Fehlerrate für echte Daten taugen sie nicht.
