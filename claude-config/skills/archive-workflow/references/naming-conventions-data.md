# Naming Conventions: Data Projects

## Overview

This reference defines naming conventions for data projects (data engineering, ETL pipelines, data warehouses, ML projects). Based on best practices from data engineering and MLOps communities.

## Data Files

### Raw Data

| Element | Convention | Example |
|---------|------------|---------|
| Location | `data/raw/{source}/` | `data/raw/api/`, `data/raw/manual/` |
| Naming | `{source}-{YYYYMMDD}-{description}.{ext}` | `salesforce-20240115-accounts.parquet` |
| Metadata | `{datafile}.meta.yaml` | `salesforce-20240115-accounts.parquet.meta.yaml` |

**Source Categories**:
- `api/` - API extracts
- `manual/` - Manual uploads
- `db/` - Database dumps
- `sftp/` - File transfers
- `streaming/` - Stream captures

### Processed Data

| Element | Convention | Example |
|---------|------------|---------|
| Location | `data/processed/` | Cleaned, transformed data |
| Naming | `{domain}-{description}-{YYYYMMDD}.{ext}` | `sales-monthly-summary-20240115.parquet` |
| Partitioned | `{domain}/{partition}/` | `sales/year=2024/month=01/` |

### Feature Stores

| Element | Convention | Example |
|---------|------------|---------|
| Location | `features/` or `data/features/` | Feature engineering outputs |
| Naming | `{entity}-{feature-group}-{version}.parquet` | `customer-demographics-v1.parquet` |
| Online | `{entity}-{feature}-latest.json` | `customer-churn-score-latest.json` |

## Notebooks

### Data Exploration

| Stage | Convention | Example |
|-------|------------|---------|
| EDA | `01-{initials}-eda-{dataset}.ipynb` | `01-da-eda-sales-data.ipynb` |
| Quality check | `02-{initials}-dq-{dataset}.ipynb` | `02-da-dq-customer-records.ipynb` |
| Profiling | `03-{initials}-profile-{dataset}.ipynb` | `03-da-profile-transactions.ipynb` |

### Model Development

| Stage | Convention | Example |
|-------|------------|---------|
| Feature engineering | `10-{initials}-features-{model}.ipynb` | `10-da-features-churn-model.ipynb` |
| Training | `20-{initials}-train-{model}.ipynb` | `20-da-train-churn-model.ipynb` |
| Evaluation | `30-{initials}-eval-{model}.ipynb` | `30-da-eval-churn-model.ipynb` |
| Inference | `40-{initials}-inference-{model}.ipynb` | `40-da-inference-churn-model.ipynb` |

## Pipeline Code

### ETL Scripts

| Element | Convention | Example |
|---------|------------|---------|
| Extract scripts | `extract_{source}.py` | `extract_salesforce.py` |
| Transform scripts | `transform_{domain}.py` | `transform_sales.py` |
| Load scripts | `load_{destination}.py` | `load_warehouse.py` |
| Full pipelines | `pipeline_{domain}.py` | `pipeline_daily_sales.py` |

### Directory Structure

```
pipelines/
├── extract/
│   ├── __init__.py
│   ├── extract_salesforce.py
│   └── extract_postgres.py
├── transform/
│   ├── __init__.py
│   ├── transform_sales.py
│   └── transform_customers.py
├── load/
│   ├── __init__.py
│   └── load_bigquery.py
└── orchestration/
    ├── dag_daily_etl.py
    └── dag_weekly_report.py
```

### Airflow DAGs

| Element | Convention | Example |
|---------|------------|---------|
| DAG files | `dag_{frequency}_{domain}.py` | `dag_daily_sales_etl.py` |
| DAG IDs | `{frequency}_{domain}_{action}` | `daily_sales_etl` |
| Task IDs | `{action}_{target}` | `extract_salesforce`, `transform_sales` |

## Models

### Model Files

| Element | Convention | Example |
|---------|------------|---------|
| Training artifacts | `{model-type}-{YYYYMMDD}-{version}.{ext}` | `xgboost-20240115-v1.pkl` |
| Serving models | `{model-name}-{version}/` | `churn-predictor-v2/` |
| Checkpoints | `checkpoint-{epoch}-{metric}.pt` | `checkpoint-050-0.95.pt` |

### Model Registry Structure

```
models/
├── churn-predictor/
│   ├── v1/
│   │   ├── model.pkl
│   │   ├── config.yaml
│   │   └── metrics.json
│   └── v2/
│       ├── model.pkl
│       ├── config.yaml
│       └── metrics.json
└── demand-forecaster/
    └── v1/
        └── ...
```

### Model Metadata

| File | Purpose | Required Fields |
|------|---------|-----------------|
| `config.yaml` | Hyperparameters | model_type, params, features |
| `metrics.json` | Performance | accuracy, f1, auc, etc. |
| `signature.yaml` | Input/output schema | input_schema, output_schema |

## Configuration Files

### Pipeline Configs

| Element | Convention | Example |
|---------|------------|---------|
| Environment | `config-{env}.yaml` | `config-prod.yaml`, `config-dev.yaml` |
| Pipeline | `{pipeline}-config.yaml` | `sales-etl-config.yaml` |
| Credentials | `{service}-credentials.yaml` | `bigquery-credentials.yaml` (gitignored) |

### Schema Files

| Element | Convention | Example |
|---------|------------|---------|
| JSON Schema | `schema-{entity}.json` | `schema-customer.json` |
| Avro | `{entity}.avsc` | `customer.avsc` |
| Protobuf | `{entity}.proto` | `customer.proto` |

## Directory Structure Template

```
project/
├── data/
│   ├── raw/
│   │   ├── api/
│   │   ├── db/
│   │   └── manual/
│   ├── processed/
│   ├── features/
│   └── external/
├── notebooks/
│   ├── exploration/
│   └── modeling/
├── pipelines/
│   ├── extract/
│   ├── transform/
│   ├── load/
│   └── orchestration/
├── models/
│   └── {model-name}/
├── configs/
│   ├── config-dev.yaml
│   └── config-prod.yaml
├── schemas/
├── tests/
├── docs/
├── README.md
├── CLAUDE.md
└── .gitignore
```

## Versioning Strategies

### Data Versioning

| Strategy | When to Use | Example |
|----------|-------------|---------|
| Date-based | Daily/periodic snapshots | `sales-20240115.parquet` |
| Partition | Large datasets | `year=2024/month=01/day=15/` |
| Hash-based | Content-addressed | DVC, LakeFS |

### Model Versioning

| Strategy | When to Use | Example |
|----------|-------------|---------|
| Semantic | Production models | `v1.2.3` |
| Date-based | Development | `20240115` |
| Git-hash | Exact reproducibility | `abc123f` |

## Anti-Patterns

| Bad | Good | Reason |
|-----|------|--------|
| `data_final_v2_REAL.csv` | `data-20240115.csv` | Use dates |
| `model.pkl` | `xgboost-20240115-v1.pkl` | Include metadata |
| `etl.py` | `extract_salesforce.py` | Be specific |
| `config.yaml` | `config-prod.yaml` | Specify environment |
| `notebook.ipynb` | `01-da-eda-sales.ipynb` | Numbered, prefixed |

## Large File Handling

### Files to Git-Track

- Code (`.py`, `.sql`)
- Configs (`.yaml`, `.json`)
- Schemas (`.avsc`, `.proto`)
- Documentation (`.md`)
- Small samples (`samples/` < 1MB each)

### Files to NOT Git-Track

- Raw data (use DVC, S3, GCS)
- Models (use model registry)
- Credentials (use secrets manager)
- Large binaries (use LFS if necessary)

## Metadata Standards

### Required Metadata for Data Files

```yaml
# data-20240115.meta.yaml
source: salesforce
extracted_at: 2024-01-15T08:30:00Z
schema_version: 1.0
row_count: 150000
columns:
  - name: customer_id
    type: string
  - name: amount
    type: float
```

### Required Metadata for Models

```yaml
# model-config.yaml
model_name: churn-predictor
version: 1.0.0
trained_at: 2024-01-15T10:00:00Z
framework: scikit-learn
features:
  - customer_tenure
  - monthly_spend
target: churn_flag
metrics:
  accuracy: 0.92
  f1_score: 0.88
```

## References

- [Cookiecutter Data Science](https://drivendata.github.io/cookiecutter-data-science/)
- [MLOps Best Practices](https://ml-ops.org/)
- [Data Engineering Cookbook](https://github.com/andkret/Cookbook)
- [DVC Documentation](https://dvc.org/doc)
