# Grafikeinheit für die Partnersuche

Auf dem MacBook Air (M5) kostet die lokale Ausrichtung der Partnersuche rund ein Viertel der
Rechenzeit (Profil vom 26.11.2026: 24 %). Es sind viele tausend voneinander unabhängige
Ausrichtungen je Block, genau die Art Arbeit, für die eine Grafikeinheit gebaut ist.
`fq2bam --grafikeinheit an` lässt sie dort rechnen. Alles bleibt auf dem Gerät: Die
Grafikeinheit sitzt im selben Rechner, beim Mac teilt sie sich den Arbeitsspeicher mit dem
Prozessor. Quelltext: `crates/gc-gpu` (wgpu: Metal auf macOS, Vulkan unter Linux, DX12 unter
Windows) und der Rechenkern `crates/gc-gpu/src/rescue.wgsl`.

## Ablauf
1. Vor dem Paaren plant jedes Paar seine Fenster ([Partnersuche gesammelt
   rechnen](read-paare.md#partnersuche-gesammelt-rechnen)).
2. Für jedes Fenster werden der Partner (in der Richtung des Fensters) und das Stück Referenz
   bereitgelegt. Fenster, deren Read länger ist, als der Byte-Lauf sicher fasst (Score über
   250), rechnet der Prozessor.
3. Die Grafikeinheit rechnet alle übrigen in einem Aufruf, eine Ausrichtung je Faden. Jeder
   Faden geht denselben gestreiften Byte-Lauf wie der Prozessor (siehe
   [Read-Paare](read-paare.md)), mit 16 Spuren als einzelnen Werten, in derselben
   Reihenfolge einschließlich der Nachkorrektur. Er liefert den besten Score, seine Zeile, die
   Read-Position und das Maximum jeder Zeile.
4. Aus den Zeilenmaxima bildet der Prozessor den zweitbesten Score (`local::from_rows`) und rechnet
   wie bisher rückwärts den Anfang des Alignments (`local::with_start`).
5. Fällt die Grafikeinheit aus, rechnet der Prozessor den Block (`GPU.FEHLER`).

## Warum die Ausgabe gleich bleibt
Der Rechenkern ist ganzzahlig und geht Schritt für Schritt den Weg des Prozessors. Ein Nachbau
dieser Schritte in Rust (`local::kernel_reference`) wird in den Tests gegen den Byte-Lauf des
Prozessors geprüft (300 zufällige Fälle mit Abweichungen, Lücken und N), ebenso ein ganzer Block
von Paaren mit diesem Nachbau als Rechenkern gegen den Block ohne. Der Rechenkern selbst wird bei
jedem Test als WGSL geprüft (Syntax und Typen). Ob die Grafikeinheit genau dasselbe liefert, prüft
`cargo run --release -p gc-gpu --example gpubench` auf dem Gerät mit Grafikeinheit: Es rechnet
zufällige Fenster auf beiden und bricht bei der ersten Abweichung mit Fehler ab.

## Messung
Steht aus: Die Cloud-Rechner der Entwicklung haben keine Grafikeinheit. Gemessen wird auf dem
MacBook Air (M5) mit `gpubench` und einem Lauf auf 1 Mio. PJ-Paaren gegen `--grafikeinheit aus`.

## Grenzen
- Bisher nur die Partnersuche. Das Verlängern (17 % auf dem Mac) folgt, wenn sich die
  Grafikeinheit hier lohnt.
- Während die Grafikeinheit rechnet, warten die Kerne des Prozessors auf das Ergebnis des Blocks.
- Reads über 256 Basen rechnet immer der Prozessor.
