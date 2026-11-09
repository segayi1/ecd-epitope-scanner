# External software: download and configuration

Use Windows for the local website and WSL2 Ubuntu 22.04 for Linux/CUDA model workers. Models, weights and scientific executables are not bundled. Obtain them from the upstream sources below and accept their licenses yourself.

Start each command in the repository root unless a different directory is shown. Replace `/mnt/e/path/to/ecd-epitope-scanner` with the WSL path to your checkout and `<linux-user>` with your own Linux username. Quote paths that contain spaces. The commands below target the Python 3.11, PyTorch 2.5.1 and CUDA 12.4 model environments used for the recorded run; the NVIDIA driver must support CUDA 12.4. Keep the following asset layout:


```text
model_assets/
  models/spurs/
    spurs/
    spurs_multi/
    esm2_t33_650M_UR50D.pt
  vendor/spurs/
  vendor/discotope3/
    models/
  vendor/foldx/
    foldx_*.exe
```

SPURS/ESM weights need about 8 GB; DiscoTope and dependencies need additional space. Allow 20–30 GB including environments and an example run. These are estimates.

## 1. WSL2 and Miniconda

In an administrator PowerShell terminal:


```powershell
wsl --install -d Ubuntu-22.04
```

Restart Windows if requested. Open Ubuntu 22.04 once and finish the initial Linux username/password setup. Then, in PowerShell:

```powershell
wsl -d Ubuntu-22.04 -- bash -lc "sudo apt update && sudo apt install -y git curl wget unzip build-essential"
```

Open Ubuntu and install Miniconda:


```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p "$HOME/miniconda3"
"$HOME/miniconda3/bin/conda" init bash
source "$HOME/.bashrc"
```

Check the Linux account and home directory. Use this account in the worker Python paths below:


```bash
whoami
echo "$HOME"
```

## 2. SPURS and weights

In WSL, clone the pinned revision, apply the supplied inference patch once, then install the worker environment and download its models:


```bash
cd "/mnt/e/path/to/ecd-epitope-scanner"
mkdir -p model_assets/vendor model_assets/models/spurs
git clone --branch fix/py311-inference-beta \
  https://github.com/luo-group/SPURS.git model_assets/vendor/spurs
# Revision used by the strict acceptance run.
git -C model_assets/vendor/spurs checkout 70cbb99af5cdf1f6da392389aaf93fd245610efb
git -C model_assets/vendor/spurs apply --check ../../../docs/patches/spurs-inference.patch
git -C model_assets/vendor/spurs apply ../../../docs/patches/spurs-inference.patch

"$HOME/miniconda3/bin/conda" create -n ecd-spurs -c conda-forge python=3.11 pip -y
"$HOME/miniconda3/bin/conda" run -n ecd-spurs python -m pip install --upgrade pip
"$HOME/miniconda3/bin/conda" run -n ecd-spurs python -m pip install \
  torch==2.5.1 --index-url https://download.pytorch.org/whl/cu124
"$HOME/miniconda3/bin/conda" run -n ecd-spurs python -m pip \
  install -c docs/constraints-spurs-py311-cu124.txt \
  -r model_assets/vendor/spurs/requirements.inference.txt
"$HOME/miniconda3/bin/conda" run -n ecd-spurs python -m pip \
  install -c docs/constraints-spurs-py311-cu124.txt -e model_assets/vendor/spurs

"$HOME/miniconda3/envs/ecd-spurs/bin/hf" download cyclization9/SPURS \
  spurs/.hydra/config.yaml spurs/checkpoints/best.ckpt \
  spurs_multi/.hydra/config.yaml spurs_multi/checkpoints/best.ckpt \
  --revision 0cc7a565af8f31eb122819f95a9d16e27b3d1596 \
  --local-dir model_assets/models/spurs
curl --fail -L \
  https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt \
  -o model_assets/models/spurs/esm2_t33_650M_UR50D.pt

# fair-esm uses the PyTorch hub checkpoint cache when loading this model.
hub_dir=$("$HOME/miniconda3/envs/ecd-spurs/bin/python" -c 'import torch; print(torch.hub.get_dir())')
mkdir -p "$hub_dir/checkpoints"
cp model_assets/models/spurs/esm2_t33_650M_UR50D.pt "$hub_dir/checkpoints/"
curl --fail -L \
  https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50D-contact-regression.pt \
  -o "$hub_dir/checkpoints/esm2_t33_650M_UR50D-contact-regression.pt"
```

The MIT patch exposes `encode_wt` / `decode_mutations` for WT representation reuse and adds the inference imports `lmdb` and `atom3d`. It has been checked against pristine code at the specified revision and reproduces the corresponding tested local files. The upstream license is retained in [patches/SPURS-LICENSE](patches/SPURS-LICENSE).

Check the CUDA environment:


```bash
"$HOME/miniconda3/envs/ecd-spurs/bin/python" -c \
  "import torch; print(torch.__version__, torch.cuda.is_available())"
```

Check both SPURS checkpoints, configuration files, the ESM weight file, and the patched API:


```bash
test -f model_assets/models/spurs/spurs/checkpoints/best.ckpt && \
test -f model_assets/models/spurs/spurs_multi/checkpoints/best.ckpt && \
test -f model_assets/models/spurs/spurs/.hydra/config.yaml && \
test -f model_assets/models/spurs/spurs_multi/.hydra/config.yaml && \
test -f model_assets/models/spurs/esm2_t33_650M_UR50D.pt && \
echo 'SPURS files present' || echo 'A required SPURS file is missing'
"$HOME/miniconda3/envs/ecd-spurs/bin/python" -c \
  "from spurs.models.stability.spurs_multi import SPURSMulti; assert callable(SPURSMulti.encode_wt); assert callable(SPURSMulti.decode_mutations); print('SPURS API OK')"
```

The ESM contact-regression checkpoint is also required in `torch.hub.get_dir()/checkpoints`, as configured above.
The pinned Hugging Face revision contains the four SPURS files used in the recorded run. Its checkpoint digests and configuration files match the recorded inputs. The constraints file pins direct inference packages observed in that run; it is not a complete lockfile for every transitive dependency.

## 3. DiscoTope3

DiscoTope3 uses CC BY-NC 4.0. Contact its authors for commercial permission. In WSL:


```bash
cd "/mnt/e/path/to/ecd-epitope-scanner"
git clone https://github.com/Magnushhoie/DiscoTope-3.0.git \
  model_assets/vendor/discotope3
git -C model_assets/vendor/discotope3 checkout 35d9f2e55f97eaba2a7acefbc394db58fb9670bc

# Align the upstream requirement with the PyTorch version used in the accepted run.
sed -i 's/^torch==2\.10\.\*$/torch==2.5.1/' model_assets/vendor/discotope3/requirements.txt
grep -Fx 'torch==2.5.1' model_assets/vendor/discotope3/requirements.txt

"$HOME/miniconda3/bin/conda" create -n ecd-discotope -c conda-forge \
  python=3.11 pip gcc gxx pybind11 -y
"$HOME/miniconda3/bin/conda" run -n ecd-discotope python -m pip install --upgrade pip
"$HOME/miniconda3/bin/conda" run -n ecd-discotope python -m pip install \
  torch==2.5.1 --index-url https://download.pytorch.org/whl/cu124
"$HOME/miniconda3/bin/conda" run -n ecd-discotope python -m pip \
  install -c docs/constraints-discotope-py311-cu124.txt \
  -r model_assets/vendor/discotope3/requirements.txt
"$HOME/miniconda3/bin/conda" run -n ecd-discotope python -m pip \
  install -c docs/constraints-discotope-py311-cu124.txt model_assets/vendor/discotope3

cd model_assets/vendor/discotope3
unzip -o models.zip
cd ../../..

# DiscoTope also loads the ESM-IF1 checkpoint from the PyTorch hub cache.
hub_dir=$("$HOME/miniconda3/envs/ecd-discotope/bin/python" -c 'import torch; print(torch.hub.get_dir())')
mkdir -p "$hub_dir/checkpoints"
curl --fail -L \
  https://dl.fbaipublicfiles.com/fair-esm/models/esm_if1_gvp4_t16_142M_UR50.pt \
  -o "$hub_dir/checkpoints/esm_if1_gvp4_t16_142M_UR50.pt"
```

The pinned upstream revision matches the accepted DiscoTope source and model archive. Its original requirements specify PyTorch 2.10; the command above updates that requirement in the local clone to match the 2.5.1 CUDA 12.4 environment used in the acceptance run. The constraints file records direct dependencies from that environment. Check the model files and ESM-IF1 cache:


```bash
test "$(find model_assets/vendor/discotope3/models -maxdepth 1 -type f -name 'XGB_*_of_100.json' | wc -l)" -eq 100 && \
test -f "$hub_dir/checkpoints/esm_if1_gvp4_t16_142M_UR50.pt" && \
echo 'DiscoTope files present' || echo 'A DiscoTope model or ESM-IF1 checkpoint is missing'
"$HOME/miniconda3/envs/ecd-discotope/bin/python" -c \
  "import torch, discotope3, torch_geometric, xgboost; print(torch.__version__, torch.cuda.is_available())"
```

## 4. Optional ProCaliper

ProCaliper provides an optional SASA/rSASA cross-check in a separate worker process. It does not determine finite-probe accessibility. Its package metadata specifies GPL-3.0-only; review upstream terms before redistribution or modification.

Install it in the main Windows virtual environment if required:


```powershell
python -m pip install procaliper==0.2.1
python -c "import procaliper; print('ProCaliper OK')"
```

## 5. Licensed FoldX 5.1

1. Register at <https://foldxsuite.crg.eu/> and accept the applicable license.
2. Download the FoldX 5.1 Windows distribution from that site.
3. Extract the supplied files into `model_assets\vendor\foldx\`, keeping the executable and its accompanying files together.
4. Set `ECD_FOLDX_BINARY` to the executable's full path in the PowerShell session that launches the application.

FoldX is obtained through the licensed download portal. Commercial use and redistribution are governed by its license.


```powershell
New-Item -ItemType Directory -Force model_assets\vendor\foldx | Out-Null
Test-Path model_assets\vendor\foldx\foldx_*.exe
```

## 6. Configure Windows worker commands

Run the following in the checkout root. The worker launcher uses the `.venv` Python created by the README installation commands; substitute your environment's Python path if you used a different location. Replace `<linux-user>` before executing it; if the Linux home directory printed in step 1 is not `/home/<linux-user>`, substitute the actual home directory in both worker Python paths. JSON argument arrays preserve paths containing spaces.


```powershell
$env:ECD_MODEL_ASSETS_ROOT = "$PWD\model_assets"
$env:ECD_FOLDX_BINARY = (Get-ChildItem "$PWD\model_assets\vendor\foldx\foldx_*.exe" |
  Select-Object -First 1).FullName

$python = (Resolve-Path '.\.venv\Scripts\python.exe').Path
$bridge = "$PWD\integrations\wsl_bridge.py"
$worker = "$PWD\integrations\workers"
$spursPython = "/home/<linux-user>/miniconda3/envs/ecd-spurs/bin/python"
$discotopePython = "/home/<linux-user>/miniconda3/envs/ecd-discotope/bin/python"

$spurs = @(
  $python, $bridge, "--distro", "Ubuntu-22.04",
  "--python", $spursPython,
  "--worker", "$worker\spurs_gpu_worker.py"
)
$env:ECD_V4_SPURS_SINGLE_COMMAND = ConvertTo-Json -InputObject $spurs -Compress
$env:ECD_V4_SPURS_MULTI_COMMAND = $env:ECD_V4_SPURS_SINGLE_COMMAND

$discotope = @(
  $python, $bridge, "--distro", "Ubuntu-22.04",
  "--python", $discotopePython,
  "--worker", "$worker\discotope_gpu_worker.py"
)
$env:ECD_V4_DISCOTOPE_COMMAND = ConvertTo-Json -InputObject $discotope -Compress

$env:ECD_V4_PROCALIPER_COMMAND = ConvertTo-Json -InputObject @(
  $python, "$worker\procaliper_worker.py"
) -Compress

$env:ECD_V4_FOLDX_COMMAND = ConvertTo-Json -InputObject @(
  $python, "$worker\foldx_worker.py"
) -Compress
```

ProCaliper is optional: if you skip its installation, omit its command setting and its unavailable cross-check will be reported. The other worker commands and asset paths must refer to the environments/files you installed.

These variables apply only to the current PowerShell session. Reapply them in a new terminal or keep a private local startup script. Do not commit local paths or credentials.

Inspect the configuration before starting the website:


```powershell
$env:ECD_MODEL_ASSETS_ROOT
$env:ECD_FOLDX_BINARY
$env:ECD_V4_SPURS_SINGLE_COMMAND
$env:ECD_V4_DISCOTOPE_COMMAND
Test-Path "$env:ECD_MODEL_ASSETS_ROOT\models\spurs\spurs\checkpoints\best.ckpt"
Test-Path "$env:ECD_FOLDX_BINARY"
```

Then launch the service from the same terminal:


```powershell
$env:ECD_DATA_DIR = "$PWD\workspace_data"
python -m ecd_epitope_scanner serve
```

Restart the service whenever configuration changes; an existing Python process does not inherit later environment changes.

## 7. Verify setup and diagnose failures

First run the README preview command to check main-package installation, mapping, and exports. Then use `strict` and a new output directory from the same configured terminal.

- Confirm `execution_protocol.json` has `name: ordinary_site_first_v2`.
- Read `recommendation_diagnostics_v41.json` for actual worker outcomes and exclusions. `tool_inventory.json` records tool identity; an installed tool may still not have been configured/invoked for a run.
- Read the per-region queues and remaining gaps in `ordinary_candidate_ledger.json`. A completed job does not imply a candidate passed every scientific gate.

| Symptom | Check |
| --- | --- |
| `TOOL_NOT_CONFIGURED` | Worker command, selected Python environment, installed package, and service restart. ProCaliper may be omitted. |
| WSL path missing | Linux username, distribution name, conda environment path, and Windows checkout accessibility from WSL. |
| ESM checkpoint missing | Both main and contact-regression weights under the SPURS environment's PyTorch hub checkpoint directory. |
| CUDA unavailable / OOM | Windows driver and `torch.cuda.is_available()` in WSL. Workers record OOM/retry/CPU fallback; CPU execution may be slower. |
| FoldX unavailable | `ECD_FOLDX_BINARY`, license status, and extracted accompanying files. |

The recorded strict acceptance used existing local model environments. The pinned revisions and direct package versions were checked against that run; these installation commands have not been repeated on a fresh machine. Record the versions, model identities, and hardware actually used for your analysis.

## 8. Interpretation

SPURS stability, DiscoTope epitope propensity, ProCaliper SASA, and FoldX ΔΔG are different forms of computational evidence. Preview is for checking workflow/geometry/exports. Strict results still require experimental validation; missing, timed-out, or invalid model evidence retains its failure status.
