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
9. [Scripts de recuperação e utilidades](#9-scripts-de-recuperação-e-utilidades)
10. [Amostra de galáxias de campo (controle)](#10-amostra-de-galáxias-de-campo-controle)
11. [Dados do S-PLUS](#11-dados-do-s-plus)
12. [Citação](#12-citação)
13. [Solução de problemas](#13-solução-de-problemas)

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

A estrutura do projeto é mostrada abaixo. `ejecutable.py` cria todos os diretórios de saída automaticamente — **só é necessário criar `Catalogos/` e colocar o catálogo de entrada antes da primeira execução**.

```
MorphoPlus/
├── Catalogos/                  # ← criar este e colocar SPLUS_Table.csv dentro
│   └── SPLUS_Table.csv         # seu catálogo de entrada
├── Field_Img/                  # criado automaticamente por ejecutable.py
│   ├── det/                    # imagens de detecção e mapas de segmentação
│   ├── mask/                   # máscaras binárias em formato FITS
│   └── psf/                    # PSFs Moffat em formato FITS por filtro
├── Out_img/                    # criado automaticamente por ejecutable.py
├── config.py                   # ← SUAS credenciais (nunca subir ao git)
├── ejecutable.py               # driver principal — gera ejecutable.sh
├── Recortar.py
├── mascara.py
├── segmetation.py
├── table_generation.py
├── leer_header_output.py
├── Img_galxgal.py
├── psf_new.py
├── count_fits.py               # verificação de progresso (galáxias prontas vs. faltantes)
├── generate_missing_fits.py    # recuperação: regera .sh com os GALFITM pendentes
├── generate_psf_from_field.py  # recuperação: reconstrói uma PSF faltante a partir de uma imagem de campo
├── gauss_5_0_9x9.conv          # Filtro de convolução para o SExtractor
├── sextopsfex.param            # Arquivo de parâmetros do SExtractor
└── galfitm-1.4.4-linux-x86_64  # Binário do GALFITM
```

```bash
mkdir Catalogos   # apenas este diretório é necessário antes de executar
```

---

## 4. Catálogo de entrada

Colocar a tabela de entrada em:
```
Catalogos/SPLUS_Table.csv
```

O catálogo deve conter as seguintes colunas:

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

### Download de dados do S-PLUS

Os catálogos fotométricos e imagens do S-PLUS estão disponíveis em **https://splus.cloud**.

---

## 5. Configuração

### Passo 0 — Configurar suas credenciais do S-PLUS

Abra `config.py` e substitua os valores de exemplo pelos seus:

```python
SPLUS_USERNAME = "seu_usuario_aqui"   # ← seu nome de usuário no S-PLUS
SPLUS_PASSWORD = "sua_senha_aqui"     # ← sua senha do S-PLUS
```

Cadastre-se em **https://splus.cloud** se ainda não tiver uma conta.

> **Segurança:** `config.py` está listado no `.gitignore` e **nunca** será enviado ao repositório. Não compartilhe este arquivo nem cole seu conteúdo em issues ou pull requests.

### Parâmetros de grade e download

Todos os parâmetros da grade são definidos **apenas em `ejecutable.py`**. Não é necessário editar nenhum outro arquivo para uma execução padrão.

| Parâmetro | Arquivo | Descrição |
|-----------|---------|-----------|
| `size` | `ejecutable.py` | Lado de cada sub-imagem em pixels. Padrão: `550`. |
| `c` | `ejecutable.py` | Lista de centros da grade (pixels). Padrão: 20 posições de 275 a 10725 em passos de 550. |
| `GALFITM_BIN` | `ejecutable.py` | Caminho/nome do binário do GALFITM. Padrão: `./galfitm-1.4.4-linux-x86_64`. |
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
- Calcula as coordenadas em pixels (X, Y) via WCS se estiverem ausentes (usando um stamp pequeno de 15×15 px na banda R — muito mais rápido que baixar o campo inteiro)
- Pré-calcula a lista explícita de comandos GALFITM iterando sobre a grade espacial
- Gera `ejecutable.sh` com todos os comandos para a execução completa

### Passo 2 — Executar o pipeline completo

```bash
chmod +wrx ejecutable.sh
./ejecutable.sh
```

O `ejecutable.sh` executa as seguintes etapas automaticamente, com logging em `morphoplus_run.log`:

1. **Etapa 1** — Baixa as 12 imagens por campo em paralelo, recorta sub-imagens na grade espacial e constrói PSFs Moffat a partir do header FITS (`Recortar.py`)
2. **Etapa 2** — Executa os scripts do SExtractor (`dopsfex_mask_*.sh`) para gerar mapas de segmentação
3. **Etapa 3** — Cria máscaras binárias e arquivos de entrada do GALFITM (`mascara.py`)
4. **Etapa 4** — Executa o GALFITM **sequencialmente** em cada sub-imagem de grupo usando uma lista de comandos explícita e pré-calculada. Falhas de ajustes individuais são logadas mas não abortam o pipeline (`set +e` está ativo nesta etapa)
5. **Etapa 5** — Executa o GALFITM galáxia por galáxia via `ejecutable_gal.sh` (apenas se existir)
6. **Etapa 6** — Lê os headers de saída do GALFITM e escreve o catálogo final de resultados (`leer_header_output.py`)

> **Por que GALFITM sequencial?** Versões anteriores do pipeline tentavam paralelizar as chamadas ao GALFITM, mas isso causava travamentos. A versão atual usa uma lista de comandos explícita e sequencial — mais lenta, mas robusta.

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
      │                                    (com lista explícita de comandos GALFITM)
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
  ├─ [Etapa 4] galfitm galfit_*_*_*.input   (sequencial, set +e)
  │     Ajuste multi-banda de Sérsic para sub-imagens de grupos
  │
  ├─ [Etapa 5] ejecutable_gal.sh   (somente se existir)
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
| `galfit_*.input` | Arquivos de entrada do GALFITM (sub-imagens de grupo) |
| `Gal_*.input` | Arquivos de entrada do GALFITM (galáxias individuais) |
| `Galfitm_*.fits` | FITS de saída do GALFITM para grupos (entrada / modelo / resíduo) |
| `Gal_*.fits` | FITS de saída do GALFITM para galáxias individuais |
| `Out_img/*.svg` | Imagens de visualização de três painéis |
| `ejecutable.sh` | Script de execução gerado automaticamente |
| `dopsfex_mask_*.sh` | Scripts do SExtractor gerados automaticamente |
| `morphoplus_run.log` | Log completo da execução (criado por `ejecutable.sh`) |

---

## 9. Scripts de recuperação e utilidades

O GALFITM roda sobre centenas de sub-imagens e galáxias individuais, e uma execução completa pode levar várias horas. Para suportar execuções longas e se recuperar de interrupções (quedas de energia, processos terminados, PSFs faltantes, ...), o pipeline inclui três utilitários pequenos.

### 9.1 `count_fits.py` — verificação de progresso

Conta quantos arquivos `.fits` de saída do GALFITM existem em disco, comparado com quantos são esperados pelo catálogo. Útil a qualquer momento para verificar o avanço do pipeline.

```bash
python count_fits.py
```

Saída típica:

```
GALFITM .fits: 312/400
  Grupos: 290/378
  Indiv:  22/22
⏳ 78.0%
```

Quando tudo está completo:

```
GALFITM .fits: 400/400
  Grupos: 378/378
  Indiv:  22/22
✅ COMPLETO
```

### 9.2 `generate_missing_fits.py` — retomar após uma queda

Se o GALFITM for interrompido no meio da execução (queda de energia, processo terminado, reinício do sistema, ...), este script verifica o diretório de trabalho buscando arquivos `.input` **sem** o correspondente `.fits` de saída e escreve um script de recuperação `missing_fits.sh` contendo **apenas** os comandos GALFITM faltantes.

```bash
python generate_missing_fits.py
chmod +x missing_fits.sh
./missing_fits.sh
```

O script de recuperação:
- Usa `set +e` para que falhas individuais não abortem a execução
- Imprime uma linha de progresso `[N/total] Completado` após cada ajuste finalizado
- Pode ser executado quantas vezes for necessário — os `.fits` já gerados são pulados

Se nada estiver faltando, o script informa e termina sem criar `missing_fits.sh`.

### 9.3 `generate_psf_from_field.py` — reconstruir uma PSF faltante

Se o `Recortar.py` não conseguiu construir uma PSF para um filtro em particular (por exemplo porque `FWHM` ou `BETA` estavam ausentes do header FITS no primeiro download, ou o arquivo foi deletado), este script procura qualquer imagem S-PLUS de campo daquele filtro em disco, lê `FWHM` e `BETA` do header e regenera a PSF Moffat faltante em `Field_Img/psf/`.

Uso:

```bash
python generate_psf_from_field.py CAMPO FILTRO
```

Exemplos:

```bash
python generate_psf_from_field.py STRIPE82-0001 F378
python generate_psf_from_field.py HYDRA-0001 R
python generate_psf_from_field.py SPLUS-n03s27 F395
```

O script:
- Procura arquivos que correspondam ao padrão `*_FILTRO.fits` no diretório de trabalho e em `Field_Img/`
- Testa vários keywords comuns do header (`FWHM`, `PSF_FWHM`, `SEEING`, ...; `BETA`, `MOFFAT_BETA`, ...)
- Escreve a nova PSF em `Field_Img/psf/psf_{CAMPO}_{FILTRO}.fits`

Após regenerar a PSF, é possível executar novamente os ajustes do GALFITM afetados com `generate_missing_fits.py`.

---

## 10. Amostra de galáxias de campo (controle)

Para processar uma **amostra de galáxias de campo** (amostra de controle, fora de aglomerados):

1. Usar `imagenes_field.py` para baixar os recortes multi-banda de cada objeto.
2. Usar `Img_galxgal.py` (função `galporgal`) para gerar os arquivos de entrada do GALFITM.

Esses scripts seguem a mesma estrutura de 12 filtros, mas trabalham sobre objetos individuais em vez de sobre a grade.

---

## 11. Dados do S-PLUS

Os catálogos fotométricos e imagens do S-PLUS estão disponíveis em **https://splus.cloud**.

---

## 12. Citação

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

## 13. Solução de problemas

| Sintoma | Causa provável | Solução |
|---------|---------------|---------|
| Pipeline interrompido no meio da execução (queda de energia, kill, reboot) | Qualquer erro | Rodar `python generate_missing_fits.py` e depois `./missing_fits.sh` para retomar **apenas** os ajustes pendentes. Ver [Seção 9.2](#92-generate_missing_fitspy--retomar-após-uma-queda). |
| Quer saber o quanto o pipeline já avançou | — | Rodar `python count_fits.py` a qualquer momento. Ver [Seção 9.1](#91-count_fitspy--verificação-de-progresso). |
| PSF não criada para um filtro | FWHM/BETA ausentes do header FITS no download inicial, ou arquivo deletado | Rodar `python generate_psf_from_field.py CAMPO FILTRO` para reconstruí-la a partir de outra imagem do mesmo filtro. Ver [Seção 9.3](#93-generate_psf_from_fieldpy--reconstruir-uma-psf-faltante). |
| `mascara.py` trava / memória esgota | Arquivos FITS não são fechados no loop | Usar `mascara_fixed.py` que emprega `with fits.open()` e `del` explícito após cada leitura de `CCDData` |
| Downloads do S-PLUS travam ou falham | Timeout de rede ou rate limit da API | Reduzir `MAX_WORKERS` em `Recortar.py` para 3–4; o pipeline tenta novamente automaticamente com backoff exponencial |
| `sex: command not found` | SExtractor não instalado ou não está no PATH | Seguir os passos de instalação da Seção 2; tentar a alternativa conda |
| Input do GALFITM tem magnitudes > 30 | Problema com flags de fotometria | O pipeline substitui automaticamente magnitudes inválidas (> 30) pelo valor da banda r |
| `KeyError: 'X'` ou `KeyError: 'Y'` | Colunas de pixels ausentes no catálogo | Rodar `ejecutable.py` primeiro; ele calcula X e Y automaticamente via WCS |
| Um ajuste do GALFITM falha mas o pipeline continua | Problema de convergência em uma sub-imagem | Esperado — a Etapa 4 roda com `set +e` para que falhas individuais não abortem a execução. Usar `count_fits.py` e `generate_missing_fits.py` para inspecionar e tentar novamente. |
