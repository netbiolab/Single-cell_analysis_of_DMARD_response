"""Fine-tune the pretrained whole-human scGPT model on historical IFN CM.

Input: paths.finetune_input_h5ad, models/whole_human/{vocab.json,args.json,best_model.pt}.
Output: outputs/finetune/{model_e*.pt,best_model.pt,run.log}.
Uses the recovered custom training loop, including masking, DDP sampling and losses.
Launch with torchrun for the recorded seven-process topology; no remote telemetry.
"""
from pathlib import Path
import os
import runpy
import sys
from pipeline.common import cli, load_config, input_path, output_dir, resolve, write_json


def main():
    parser = cli(__doc__)
    parser.add_argument('--write-config-only', action='store_true')
    args = parser.parse_args()
    cfg = load_config(args.config)
    out = output_dir(cfg, 'finetune')
    rank = int(os.environ.get('RANK', 0))
    h = dict(cfg['finetune']['hparams'])
    h['load_model'] = str(resolve(cfg, cfg['paths']['pretrained_model_dir']))
    data = (output_dir(cfg, 'prepared') / 'IFN_CM_all_genes.h5ad' if cfg['finetune']['use_prepared_input'] else input_path(cfg, 'finetune_input_h5ad'))
    run = {'save_dir':str(out), 'train_adata_path':str(data), 'batch_id':'batch.id', 'celltype_id':'celltype.id', 'data_is_raw':False, 'wandb_project':'local-only'}
    # Per-rank files avoid parallel writers modifying the same JSON.
    run_path, hp_path = out / f'run_config_rank{rank}.json', out / f'hparams_rank{rank}.json'
    write_json(run_path, run)
    write_json(hp_path, h)
    if args.write_config_only:
        return
    sys.argv = ['finetune', '--config', str(run_path), '--hparams', str(hp_path)]
    runpy.run_path(str(Path(__file__).parent / 'pipeline/finetune.py'), run_name='__main__')


if __name__ == '__main__':
    main()
