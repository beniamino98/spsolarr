<!-- pd-meta: {"date-update":"2026-10-29 03:56:53"} -->
# spsolarr

`spsolarr` estende [`solarr`](https://github.com/beniamino98/solarr) con modelli
spaziali per la radiazione solare, strutture di correlazione e metodi di
interpolazione. Il repository pubblico ufficiale è
[solarr-project-public/packages/spsolarr](https://gitlab.com/solarr-project-public/packages/spsolarr).

## Autore

Beniamino Sartini - beniamino.sartini2@unibo.it.

## Requisiti e installazione

È richiesto R 4.4.0 o successivo. Le dipendenze sono dichiarate in
[DESCRIPTION](DESCRIPTION); `solarr` deve essere installato dal suo repository
ufficiale prima di installare `spsolarr`.

```r
install.packages("remotes")
remotes::install_github("beniamino98/solarr", upgrade = "never")
remotes::install_gitlab(
  "solarr-project-public/packages/spsolarr",
  upgrade = "never",
  build_vignettes = FALSE
)
```

Il codice di installazione identifica i repository ufficiali; il mantenimento
delle versioni delle dipendenze è a carico dell'ambiente R dell'utente.

## Primo utilizzo

```r
library(spsolarr)
havDistance(43.3, 12.1, 43.4, 12.2)
help(package = "spsolarr")
```

L'esempio calcola una distanza in chilometri fra due punti. Le pagine di aiuto
includono i modelli di correlazione, il clear sky spaziale e le varianti di
kriging esposte dal pacchetto.

## Struttura

- `R/`: implementazioni e funzioni esportate.
- `man/`: documentazione R.
- `tests/`: test disponibili nel progetto.
- `vignettes/`: notebook e materiali di esempio.

## Licenza e segnalazioni

Il progetto usa la licenza GPL-3, come dichiarato in `DESCRIPTION`.
Le segnalazioni possono essere aperte nelle
[issue del repository pubblico](https://gitlab.com/solarr-project-public/packages/spsolarr/-/issues).
