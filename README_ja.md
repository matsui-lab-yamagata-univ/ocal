# ocal: 有機半導体分子の軌道テール割合（OTF）計算プログラム
[![Python](https://img.shields.io/badge/python-3.11%20or%20newer-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![PyPI](https://img.shields.io/pypi/v/yu-ocal)](https://pypi.org/project/yu-ocal/)

[English](https://github.com/matsui-lab-yamagata-univ/ocal/blob/main/README.md) / 日本語

# 概要
`ocal` は有機半導体分子 1 分子の軌道テール割合（orbital tail fraction, OTF）を計算するツールです。分子構造または結晶構造（CIF）から Gaussian の入力を生成し、構造最適化を実行し、チェックポイントファイルを変換して、フロンティア軌道 4 本（NHOMO, HOMO, LUMO, NLUMO）の分子軌道 cube ファイルを生成します。各軌道について、分子のファンデルワールス表面より**外側**にあるボクセルの確率密度を積分し、その割合を OTF として報告します。

OTF は、フロンティア軌道がどれだけ分子のファンデルワールス体積からはみ出しているか、すなわち分子間の軌道の重なりにどれだけ寄与しうるかを定量化した指標です。

<p align="center">
  <img src="https://raw.githubusercontent.com/matsui-lab-yamagata-univ/ocal/main/assets/OTF.gif" alt="軌道テール割合（OTF）の概念図" width="400">
</p>

# 必要環境
* Python 3.11 以降
* NumPy 2.0.2 以降
* Pandas 2.3.3 以降

## 量子化学計算ツール
以下が必要です:
* Gaussian 16（`g16`, `formchk`, `cubegen` が `$PATH` から実行できること）

# 注意事項
* Gaussian 実行ファイルのパスを設定する必要があります。`ocal` は `g16`, `formchk`, `cubegen` をサブプロセスとして呼び出します。
* 分子は常に**中性閉殻一重項**として扱われます。電荷とスピン多重度は `0 1` に固定され、入力ファイルに含まれる電荷・スピン情報は無視されます。
* 汎関数と基底関数系は、過去に報告された OTF の値との比較可能性を保つため **B3LYP/6-31G(d,p)** に固定されています。コマンドラインオプションでは変更できません。
* 構造ファイルを入力とした場合、cube ファイルの生成前に必ず `Opt=Tight` で構造最適化が実行されます。

# インストール
`ocal` は PyPI では **`yu-ocal`** という名前で配布されています（import 名・コマンド名はいずれも `ocal` です）。

```bash
pip install yu-ocal
```

NumPy と Pandas は自動的にインストールされます。Gaussian 16 は含まれないため、別途インストールが必要です。

## インストールの確認

インストール後、以下を実行して確認できます:

```bash
ocal --help
```

# ocal 使用マニュアル

## 基本的な使い方

```bash
ocal <filename> [options]
```

### 必須引数

- `file`: 入力ファイルのパス。

`ocal` は 3 種類の入力を受け付け、それぞれに対応する段階からパイプラインを開始します。

| 入力 | 拡張子 | 実行される処理 |
|------|--------|----------------|
| 構造ファイル | `.gjf`, `.com`, `.xyz`, `.mol`, `.mol2`, `.cif` | gjf 生成 → Gaussian（Opt=Tight）→ formchk → cubegen → OTF |
| フォーマット済みチェックポイント | `.fchk` | cubegen → OTF |
| cube ファイル | `.cube` | OTF のみ |

`.fchk` および `.cube` を入力とする場合は `-s, --skip-gaussian` が**必須**です。逆に、構造ファイルを入力とする場合は指定できません。

> **注意:** `.cif` を入力とした場合、使用されるのは**最初の独立分子**（インデックス 0）のみです。単位格子内に複数の独立分子が存在する場合（`Z' > 1`）は、その旨のメッセージが表示されます。

### 基本的な例

```bash
# 分子構造ファイルから一括で計算
ocal xxx.xyz

# 既存の Gaussian フォーマット済みチェックポイントファイルから開始
ocal xxx.fchk -s

# 既に存在する cube ファイル 1 本から OTF を計算
ocal xxx_HOMO.cube -s
```

## オプション

### 計算設定

#### `-c, --cpu <number>`
Gaussian および `cubegen` が使用する CPU 数を指定します。
- **デフォルト**: `4`
- **例**: `ocal xxx.xyz -c 8`

#### `-m, --mem <memory>`
Gaussian が使用するメモリ量を GB 単位で指定します。
- **デフォルト**: `10`
- **例**: `ocal xxx.xyz -m 16`

### 計算の制御

#### `-s, --skip-gaussian`
Gaussian 計算をスキップし、既存の結果を再利用します。入力が `.fchk` または `.cube` の場合は必須で、それ以外の入力では指定できません。
- **デフォルト**: 無効
- **例**: `ocal xxx.fchk -s`

## 実用的な使用例

### 基本的な計算
```bash
# デフォルト実行（Gaussian 16、B3LYP/6-31G(d,p)、Opt=Tight）
ocal xxx.xyz

# 8 CPU・16 GB メモリを使用
ocal xxx.mol -c 8 -m 16

# 結晶構造から最初の独立分子を取り出して計算
ocal xxx.cif
```

### 結果の再利用
```bash
# 既存の fchk から cube を再生成して OTF を計算し直す
ocal xxx.fchk -s

# 既存の cube ファイル 1 本から OTF を計算し直す
ocal xxx_LUMO.cube -s
```

## 出力

### 標準出力
`ocal` は入力ファイル名、タイムスタンプ、実行した外部コマンドを表示したのち、結果を表形式で出力します:

```
----------------------------------------
orbital               OTF    density_sum
----------------------------------------
NHOMO        0.1043821735   0.9998672314
HOMO         0.1187456210   0.9998913057
LUMO         0.1352907441   0.9998745092
NLUMO        0.1490233866   0.9998501773
----------------------------------------
```

- **`OTF`**: 軌道テール割合。各原子のファンデルワールス半径より外側に位置するボクセルについて、規格化された確率密度を合計した値です。
- **`density_sum`**: cube グリッド全体で積分した規格化密度の合計。`1.0` に近い値になるはずで、1 から大きく外れている場合は cube グリッドが軌道を収めるには小さすぎるか粗すぎることを意味し、その cube から得られた OTF は信頼できません。

入力が cube ファイル 1 本の場合、行のラベルは軌道名ではなく cube ファイルのステム名になります。

### 生成されるファイル
生成されるファイルはすべて入力ファイルと同じディレクトリに、同じベース名で出力されます:

```
<入力ディレクトリ>/
├── <NAME>.gjf          # 生成された Gaussian 入力（Opt=Tight, B3LYP/6-31G(d,p)）
├── <NAME>.log          # Gaussian 出力（Windows 版では .out）
├── <NAME>.chk          # Gaussian チェックポイントファイル
├── <NAME>.fchk         # フォーマット済みチェックポイントファイル（formchk）
├── <NAME>_NHOMO.cube   # 分子軌道 cube ファイル（cubegen）
├── <NAME>_HOMO.cube
├── <NAME>_LUMO.cube
└── <NAME>_NLUMO.cube
```

#### 出力ファイルの命名
`.gjf` / `.com` を入力とした場合、入力ファイルを上書きしないようにベース名は `<NAME>_ocal` になります。cube ファイルはベース名に軌道名を付した `<NAME>_NHOMO`, `<NAME>_HOMO`, `<NAME>_LUMO`, `<NAME>_NLUMO` という名前です。

`cubegen` に渡される MO のインデックスは、fchk に記録された alpha 電子数 `na` から決定されます（NHOMO = `na-1`, HOMO = `na`, LUMO = `na+1`, NLUMO = `na+2`）。

## 補足

1. **計算時間**: 実行時間のほとんどは Gaussian の構造最適化であり、原子数の増加に伴って急激に増大します。
2. **メモリ使用量**: 大きな分子には十分なメモリを確保してください（`-m`）。
3. **Gaussian のインストール**: Gaussian 16 が必要です。Gaussian に同梱される `formchk` と `cubegen` も実行可能である必要があります。
4. **cube グリッド**: cube ファイルは `cubegen` の `-2`（fine）グリッドとヘッダーオプション `h` で生成されます。OTF の値を使う前に必ず `density_sum` を確認してください。

## トラブルシューティング

### Gaussian が正常終了しない場合
Gaussian のログに正常終了の行が無い場合、`ocal` は処理を中断します。入力ファイルと同じ場所にある `.log`（Windows 版では `.out`）を確認し、SCF や構造最適化の問題を解決してから再実行してください。

### `Failed to execute g16 / formchk / cubegen` と表示される場合
実行ファイルが見つかっていません。Gaussian の環境設定スクリプトを読み込むなどして、`g16`, `formchk`, `cubegen` が `$PATH` に含まれる状態にしてください。

### `density_sum` が 1.0 から大きく外れる場合
cube グリッドが軌道全体を捉えられていません。より広い、あるいはより細かいグリッドで cube ファイルを再生成し、その結果に対して `-s` を付けて再実行すれば、Gaussian を再計算せずに OTF を求められます。

```bash
ocal xxx.cube -s
```

### 構造ファイルが読み込めない場合
構造ファイルには様々な形式があり、ocal で読み込めないものもあります。以下をお試しください:

1. **別のソフトウェアで形式を変換する**: [Mercury](https://www.ccdc.cam.ac.uk/solutions/software/mercury/) や Open Babel などのソフトウェアでファイルを開き、再エクスポートすると解決する場合があります。
2. **お問い合わせ**: 読み込めないファイルを下記のメールアドレスにお送りいただければ、対応を検討します。

> **注意:** `.gjf` / `.com` 入力では、`元素記号 x y z` の直交座標形式のみをサポートしています。元素記号の代わりに原子番号を用いた記法、freeze フラグ、Z-matrix、ONIOM 層は受け付けません。

# 著者
[松井研究室, 有機エレクトロニクス研究センター（ROEL）, 山形大学](https://matsui-lab.yz.yamagata-u.ac.jp/)  
岡田 智悠, 尾沢 昂輝, 松井 弘之  
Email: h-matsui[at]yz.yamagata-u.ac.jp  
[at] を @ に置き換えてください
