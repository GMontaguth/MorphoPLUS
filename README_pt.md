# MorphoPlus

🌐 [English](README.md) · [Español](README_es.md) · **Português**

---

**MorphoPlus** é um pipeline semi-automático para ajuste morfológico multi-banda de galáxias em campos de aglomerados do S-PLUS usando o [GALFITM](https://www.nottingham.ac.uk/astronomy/megamorph/). O programa baixa imagens do S-PLUS, constrói PSFs, cria máscaras de segmentação e gera arquivos de entrada para o GALFITM — tudo organizado em torno de uma grade espacial que subdivide cada campo do S-PLUS em sub-imagens gerenciáveis.

---

## Índice

1. [Dependências](#1-dependências)
2. [Instalação](#2-instalação)
3. [Estrutura de diretórios](#3-estrutura-de-diretórios)
4. [Catálogo de entrada](#4-catálogo-de-entrada)
5. [Configuração](#5-configuração)
6. [Execução do pipeline](#6-execução-do-pipeline)
7. [Descrição do pipeline](#7-descrição-do-pipeline)
8. [Arquivos de saída](#8-arquivos-de-saída)
9. [Amostra de galáxias de campo (controle)](#9-amostra-de-galáxias-de-campo-controle)
10. [Query SQL do S-PLUS](#10-query-sql-do-s-plus)
11. [Citação](#11-citação)
12. [Solução de problemas](#12-solução-de-problemas)

---

## 1. Dependências

### Pacotes Python

```bash
pip install astropy numpy matplotlib splusdata
```

### Ferramentas externas

| Ferramenta | Versão | Função |
|------------|--------|--------|
| [SExtractor](https://github.com/astromatic/sextractor) | ≥ 2.25 | Detecção de fontes e mapas de segmentação |
| [PSFEx](https://github.com/astromatic/psfex) | ≥ 3.21 | Modelagem da PSF |
| [GALFITM](https://www.nottingham.ac.uk/astronomy/megamorph/) | 1.4.4 | Ajuste multi-banda de perfis de Sérsic |

---

## 2. Instalação

### Passo 1 — Bibliotecas do sistema

```bash
sudo apt-get update
sudo apt-get install libatlas-base-dev libblas-dev liblapack-dev   # ATLAS 3.6
sudo apt-get install libfftw3-dev                                  # FFTW 3.3.8
sudo apt-get install libplplot-dev                                 # PLPlot 5.9
```

### Passo 2 — SExtractor

```bash
git clone https://github.com/astromatic/sextractor.git
cd sextractor
sh autogen.sh
./configure
make -j
sudo make install
```

Verificar a instalação:
```bash
sex --version
```

> **Alternativa (conda):**
> ```bash
> conda install -c conda-forge astromatic-source-extractor
> ```
> Mais informações: https://anaconda.org/conda-forge/astromatic-source-extractor

### Passo 3 — PSFEx

```bash
git clone https://github.com/astromatic/psfex.git
cd psfex
sh autogen.sh
./configure
make -j
sudo make install
```

Verificar a instalação:
```bash
psfex --version
```

### Passo 4 — GALFITM

Baixar o binário para sua plataforma em:
https://www.nottingham.ac.uk/astronomy/megamorph/

Colocar o binário (`galfitm-1.4.4-linux-x86_64` ou equivalente) no **diretório raiz do MorphoPlus**.

---

## 3. Estrutura de diretórios

Antes de executar o pipeline, criar os seguintes diretórios vazios:

```
MorphoPlus/
├── Catalogos/                  # Catálogo de entrada e tabelas de saída
├── Field_Img/                  # Imagens baixadas e processadas
│   ├── det/                    # Imagens de detecção e mapas de segmentação
│   ├── mask/                   # Máscaras binárias em formato FITS
│   └── psf/                    # PSFs Moffat em formato FITS por filtro
├── Out_img/                    # Imagens de visualização de saída (SVG)
├── ejecutable.py
├── Recortar.py
├── mascara.py
├── segmetation.py
├── table_generation.py
├── leer_header_output.py
├── Img_galxgal.py
├── psf_new.py
├── gauss_5_0_9x9.conv          # Filtro de convolução para o SExtractor
├── sextopsfex.param            # Arquivo de parâmetros do SExtractor
└── galfitm-1.4.4-linux-x86_64 # Binário do GALFITM
```

Comando rápido de configuração:
```bash
mkdir -p Catalogos Field_Img/det Field_Img/mask Field_Img/psf Out_img
```

---

## 4. Catálogo de entrada

Colocar a tabela de entrada em:
```
Catalogos/SPLUS_Table.csv
```

O catálogo deve estar **corrigido pela extinção** e conter as seguintes colunas:

| Coluna | Descrição |
|--------|-----------|
| `ID` | Identificador único da fonte |
| `Field` | Nome do campo S-PLUS (ex. `SPLUS-S14-S01`) |
| `ra` / `dec` | Ascensão reta e declinação (graus) |
| `RA` / `DEC` | Igual ao anterior (usado para recortes) |
| `X` / `Y` | Coordenadas em pixels dentro do campo S-PLUS |
| `ELONGATION` | Razão entre os semi-eixos maior e menor (A/B) |
| `THETA` | Ângulo de posição do SExtractor (graus) |
| `FLUX_RADIUS_50` | Raio de meia luz em pixels |
| `FLUX_RADIUS_90` | Raio que contém 90% da luz em pixels |
| `r_auto` | Magnitude AUTO na banda r |
| `u_auto`, `g_auto`, `i_auto`, `z_auto` | Magnitudes AUTO nos filtros de banda larga |
| `J0378_auto` … `J0861_auto` | Magnitudes AUTO nos filtros estreitos do S-PLUS |

### Geração do catálogo com o notebook Jupyter

O notebook `splus-table-morfoplus.ipynb` automatiza a preparação do catálogo a partir de um arquivo parquet ou CSV bruto do S-PLUS. Ele renomeia colunas, calcula `ELONGATION` a partir de A/B e gera o `SPLUS_Table.csv` no formato correto.

Se o catálogo não contiver as colunas `X`, `Y` ou `ID`, o `ejecutable.py` as calculará automaticamente usando o WCS dos frames baixados.

### Query SQL para baixar fotometria do S-PLUS (DR4/DR5)

A query abaixo recupera as colunas morfológicas e fotométricas necessárias. Substituir `ra`, `dec` e `r_g` pelas coordenadas do centro do aglomerado e pelo raio de busca (em graus):

```sql
SELECT
  det.ID, det.Field, det.RA, det.DEC, det.X, det.Y,
  det.ELONGATION, det.ELLIPTICITY, det.THETA, det.A, det.B,
  det.FLUX_RADIUS_50, det.FLUX_RADIUS_90,
  u.u_petro, u.e_u_petro,
  J0378.J0378_petro, J0378.e_J0378_petro,
  J0395.J0395_petro, J0395.e_J0395_petro,
  J0410.J0410_petro, J0410.e_J0410_petro,
  J0430.J0430_petro, J0430.e_J0430_petro,
  g.g_petro, g.e_g_petro,
  J0515.J0515_petro, J0515.e_J0515_petro,
  r.r_petro, r.e_r_petro,
  J0660.J0660_petro, J0660.e_J0660_petro,
  i.i_petro, i.e_i_petro,
  J0861.J0861_petro, J0861.e_J0861_petro,
  z.z_petro, z.e_z_petro,
  pz.zml, pz.odds
FROM idr4_dual.idr4_detection_image AS det
  JOIN idr4_dual.idr4_dual_u      AS u      ON u.ID      = det.ID
  JOIN idr4_dual.idr4_dual_J0378  AS J0378  ON J0378.ID  = det.ID
  JOIN idr4_dual.idr4_dual_J0395  AS J0395  ON J0395.ID  = det.ID
  JOIN idr4_dual.idr4_dual_J0410  AS J0410  ON J0410.ID  = det.ID
  JOIN idr4_dual.idr4_dual_J0430  AS J0430  ON J0430.ID  = det.ID
  JOIN idr4_dual.idr4_dual_g      AS g      ON g.ID      = det.ID
  JOIN idr4_dual.idr4_dual_J0515  AS J0515  ON J0515.ID  = det.ID
  JOIN idr4_dual.idr4_dual_r      AS r      ON r.ID      = det.ID
  JOIN idr4_dual.idr4_dual_J0660  AS J0660  ON J0660.ID  = det.ID
  JOIN idr4_dual.idr4_dual_i      AS i      ON i.ID      = det.ID
  JOIN idr4_dual.idr4_dual_J0861  AS J0861  ON J0861.ID  = det.ID
  JOIN idr4_dual.idr4_dual_z      AS z      ON z.ID      = det.ID
  JOIN idr4_vacs.idr4_photoz      AS pz     ON pz.ID     = det.ID
WHERE 1 = CONTAINS(
  POINT('ICRS', det.RA, det.DEC),
  CIRCLE('ICRS', ra, dec, r_g)
)
```

---

## 5. Configuração

Todos os parâmetros da grade são definidos **apenas em `ejecutable.py`**. Não é necessário editar nenhum outro arquivo para uma execução padrão.

| Parâmetro | Arquivo | Descrição |
|-----------|---------|-----------|
| `size` | `ejecutable.py` | Lado de cada sub-imagem em pixels. Padrão: `550`. |
| `c` | `ejecutable.py` | Lista de centros da grade (pixels). Padrão: 20 posições de 275 a 10725 em passos de 550. |
| `MAX_WORKERS` | `Recortar.py` | Threads paralelas para download dos filtros. Padrão: `6`. Reduzir para 3–4 se a API do S-PLUS retornar erros de rate limit. |
| `DR` | `Recortar.py` | Data release preferido do S-PLUS. Padrão: `"dr6"`. Fallbacks são testados automaticamente. |
| `MAX_RETRIES` | `Recortar.py` | Tentativas máximas de download por filtro. Padrão: `5`. |

---

## 6. Execução do pipeline

### Passo 1 — Gerar o script de execução

```bash
python ejecutable.py
```

Este script:
- Lê `Catalogos/SPLUS_Table.csv`
- Calcula as coordenadas em pixels (X, Y) via WCS se estiverem ausentes
- Gera `ejecutable.sh` com todos os comandos para a execução completa

### Passo 2 — Executar o pipeline completo

```bash
./ejecutable.sh
```

O `ejecutable.sh` executa a seguinte sequência automaticamente:

1. Baixa as 12 imagens por campo em paralelo
2. Recorta sub-imagens na grade espacial
3. Constrói PSFs Moffat a partir dos keywords do header FITS
4. Executa o SExtractor para gerar mapas de segmentação
5. Cria máscaras binárias (`mascara.py`)
6. Gera os arquivos de entrada do GALFITM (`.input`)
7. Executa o GALFITM em cada sub-imagem

### Passo 3 — Ler os resultados do GALFITM

Este passo é executado automaticamente no final do `ejecutable.sh`, mas também pode ser rodado de forma independente:

```bash
python leer_header_output.py
```

Gera:
- `Catalogos/GalfitM_output.csv` — parâmetros ajustados (posição, Re, n, b/a, AP, magnitude) para as 12 bandas com incertezas
- `Out_img/*.svg` — imagens de três painéis (entrada / modelo / resíduo)

---

## 7. Descrição do pipeline

```
SPLUS_Table.csv
      │
      ▼
ejecutable.py ──────────────────────────── gera ejecutable.sh
      │
      ▼
./ejecutable.sh
  │
  ├─ [Etapa 1] Recortar.py
  │     Baixa 12 filtros (paralelo), recorta grade, constrói PSFs
  │     Gera: dopsfex_mask_{campo}.sh
  │
  ├─ [Etapa 2] dopsfex_mask_*.sh
  │     Executa SExtractor → mapas de segmentação (.seg.fits)
  │
  ├─ [Etapa 3] mascara.py
  │     Máscaras binárias + inputs do GALFITM + ejecutable_gal.sh
  │
  ├─ [Etapa 4] galfitm galfit_*_*_*.input
  │     Ajuste multi-banda de Sérsic para sub-imagens de grupos
  │
  ├─ [Etapa 5] ejecutable_gal.sh  (somente se existir)
  │     GALFITM galáxia por galáxia
  │
  └─ [Etapa 6] leer_header_output.py
        GalfitM_output.csv + imagens SVG
```

### Grade espacial

Como o GALFITM não consegue ajustar um campo inteiro do S-PLUS (≈ 10800 × 10800 px), o campo é dividido em sub-imagens centradas em uma grade regular. A configuração padrão usa 20 posições por eixo (400 sub-imagens por campo), cada uma de 550 × 550 px. Galáxias dentro de 25 px da borda de uma sub-imagem são excluídas daquele ajuste para evitar efeitos de borda.

### Convenção de máscaras

| Valor do pixel | Significado |
|----------------|-------------|
| `0` | Galáxia do grupo — **não** mascarada, o GALFITM irá ajustá-la |
| `1` | Outra fonte — mascarada, excluída do ajuste |

### Pontos-zero fotométricos

Os pontos-zero são recuperados espacialmente da API do S-PLUS (`conn.get_zp`) para cada campo, banda e posição no céu. Usa DR6 por padrão, que fornece correções de ZP por posição levando em conta os gradientes de iluminação do CCD.

### Modo galáxia individual

Para galáxias na lista `Catalogos/g_S.csv` (detectadas como isoladas ou em grupos dispersos), o `Img_galxgal.py` baixa stamps de 200 × 200 px diretamente via `conn.stamp()` e ajusta cada galáxia individualmente em vez de usar a grade.

---

## 8. Arquivos de saída

| Arquivo | Descrição |
|---------|-----------|
| `Catalogos/SPLUS_Table.csv` | Catálogo de entrada (atualizado com X, Y, ID se estavam ausentes) |
| `Catalogos/g_S.csv` | Sub-catálogo de galáxias processadas individualmente |
| `Catalogos/procesados.txt` | Log de campos já processados (para retomar execuções) |
| `Catalogos/procesados_gal.txt` | Log de galáxias individuais já processadas |
| `Catalogos/GalfitM_output.csv` | Parâmetros ajustados pelo GALFITM para todas as galáxias e bandas |
| `Field_Img/*.fits` | Recortes de sub-imagens (um arquivo por posição e filtro) |
| `Field_Img/det/*.fits` | Imagens de detecção (soma G+R+Z) e mapas de segmentação |
| `Field_Img/mask/*.fits` | Máscaras binárias |
| `Field_Img/psf/*.fits` | PSFs Moffat |
| `galfit_*.input` | Arquivos de entrada do GALFITM |
| `Galfitm_*.fits` | FITS de saída do GALFITM (entrada / modelo / resíduo) |
| `Out_img/*.svg` | Imagens de visualização de três painéis |
| `ejecutable.sh` | Script de execução gerado automaticamente |
| `dopsfex_mask_*.sh` | Scripts do SExtractor gerados automaticamente |

---

## 9. Amostra de galáxias de campo (controle)

Para processar uma **amostra de galáxias de campo** (amostra de controle, fora de aglomerados):

1. Usar `imagenes_field.py` para baixar os recortes multi-banda de cada objeto.
2. Usar `Img_galxgal.py` (função `galporgal`) para gerar os arquivos de entrada do GALFITM.

Esses scripts seguem a mesma estrutura de 12 filtros, mas trabalham sobre objetos individuais em vez de sobre a grade.

---

## 10. Query SQL do S-PLUS

Ver a [Seção 4](#4-catálogo-de-entrada) para a query SQL completa para baixar catálogos fotométricos do banco de dados do S-PLUS (DR4/DR5). A query retorna magnitudes petrosianas nas 12 bandas junto com parâmetros morfológicos e redshifts fotométricos.

---


## 11. Citação

Se você usar o MorphoPlus em sua pesquisa, por favor cite o seguinte artigo:

> Montaguth, G. P., O'Mill, A. L., Mendes de Oliveira, C., et al. (2026).
> *Galaxy Evolution in Compact Groups. III. Structural Analysis of Galaxies and Dynamical State of Non-isolated Compact Groups.*
> The Astrophysical Journal, 998(1), 91.
> https://doi.org/10.3847/1538-4357/ae2bdd

Entrada BibTeX:

```bibtex
@article{Montaguth_2026,
  doi       = {10.3847/1538-4357/ae2bdd},
  url       = {https://doi.org/10.3847/1538-4357/ae2bdd},
  year      = {2026},
  month     = {feb},
  publisher = {The American Astronomical Society},
  volume    = {998},
  number    = {1},
  pages     = {91},
  author    = {Montaguth, Gissel P. and O'Mill, Ana Laura and
               Mendes de Oliveira, Claudia and Lima-Dias, Ciria and
               Torres-Flores, Sergio and Monachesi, Antonela and
               Olave-Rojas, D. E. and Pallero, Diego and
               Humire, Pedro K. and Demarco, Ricardo and
               Telles, Eduardo and Lopes, Paulo A. A. and
               Panda, Swayamtrupta and Haack, Rodrigo F. and
               Lopes, Amanda R. and Alvarez-Candal, Alvaro and
               Smith Castelli, Analia V. and Kanaan, Antonio and
               Ribeiro, Tiago and Schoenell, William},
  title     = {Galaxy Evolution in Compact Groups. {III}. Structural Analysis
               of Galaxies and Dynamical State of Non-isolated Compact Groups},
  journal   = {The Astrophysical Journal}
}
```

---

## 12. Solução de problemas

| Sintoma | Causa provável | Solução |
|---------|---------------|---------|
| `mascara.py` trava / memória esgota | Arquivos FITS não são fechados no loop | Usar `mascara_fixed.py` que emprega `with fits.open()` e `del` explícito após cada leitura de `CCDData` |
| Downloads do S-PLUS travam ou falham | Timeout de rede ou rate limit da API | Reduzir `MAX_WORKERS` em `Recortar.py` para 3–4; o pipeline tenta novamente automaticamente com backoff exponencial |
| `sex: command not found` | SExtractor não instalado ou não está no PATH | Seguir os passos de instalação da Seção 2; tentar a alternativa conda |
| Input do GALFITM tem magnitudes > 30 | Problema com flags de fotometria | O pipeline substitui automaticamente magnitudes inválidas (> 30) pelo valor da banda r |
| Pipeline interrompido no meio da execução | Qualquer erro | Rodar novamente `python ejecutable.py` e `./ejecutable.sh`; campos e galáxias já processados são pulados via os logs em `Catalogos/` |
| `KeyError: 'X'` ou `KeyError: 'Y'` | Colunas de pixels ausentes no catálogo | Rodar `ejecutable.py` primeiro; ele calcula X e Y automaticamente via WCS |
| PSF não criada para um filtro | FWHM/BETA não encontrados no header FITS | Um aviso é impresso; o GALFITM rodará sem essa PSF — verificar a versão do DR |
