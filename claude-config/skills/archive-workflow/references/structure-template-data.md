# Structure Template: Data Projects

## Overview

This reference defines directory structure templates for data projects (data engineering, ETL pipelines, data platforms, analytics). Based on data engineering best practices and modern data stack patterns.

## Standard Data Project

### Full Template

```
project/
├── data/
│   ├── raw/                  # Unprocessed source data
│   │   ├── {source}/         # Organized by source
│   │   └── README.md
│   ├── processed/            # Cleaned, transformed data
│   │   └── {domain}/
│   ├── features/             # Feature store (if ML)
│   └── external/             # Third-party data
├── notebooks/
│   ├── exploration/          # EDA notebooks
│   │   └── 01-{initials}-eda-{dataset}.ipynb
│   └── analysis/             # Analysis notebooks
│       └── 01-{initials}-{analysis}.ipynb
├── pipelines/                # ETL/ELT code
│   ├── extract/
│   │   ├── __init__.py
│   │   └── extract_{source}.py
│   ├── transform/
│   │   ├── __init__.py
│   │   └── transform_{domain}.py
│   └── load/
│       ├── __init__.py
│       └── load_{destination}.py
├── dags/                     # Orchestration (Airflow)
│   └── dag_{frequency}_{domain}.py
├── models/                   # dbt models (if using dbt)
│   ├── staging/
│   ├── intermediate/
│   └── marts/
├── configs/
│   ├── config-dev.yaml
│   ├── config-prod.yaml
│   └── connections.yaml      # Connection configs (gitignored secrets)
├── schemas/                  # Data schemas
│   ├── raw/
│   └── processed/
├── tests/
│   ├── unit/
│   ├── integration/
│   └── data/                 # Data quality tests
├── docs/
│   ├── architecture.md
│   └── data-catalog/
├── scripts/
│   ├── setup.sh
│   └── seed.sh
├── README.md
├── CLAUDE.md
├── pyproject.toml
└── .gitignore
```

### Minimal Structure

```
project/
├── data/
│   ├── raw/
│   └── processed/
├── notebooks/
├── pipelines/
│   └── etl.py
├── configs/
│   └── config.yaml
├── README.md
├── CLAUDE.md
└── .gitignore
```

## dbt Project Structure

### Standard dbt Layout

```
project/
├── analyses/                 # Ad-hoc analyses
├── dbt_packages/             # dbt packages (gitignored)
├── logs/                     # dbt logs (gitignored)
├── macros/                   # dbt macros
│   └── generate_schema_name.sql
├── models/
│   ├── staging/              # Source cleaning
│   │   ├── {source}/
│   │   │   ├── __{source}__models.yml
│   │   │   ├── __{source}__sources.yml
│   │   │   └── stg_{source}__{table}.sql
│   │   └── README.md
│   ├── intermediate/         # Business logic
│   │   └── int_{domain}__{description}.sql
│   └── marts/                # Final tables
│       ├── core/
│       │   └── dim_{entity}.sql
│       │   └── fct_{event}.sql
│       └── marketing/
│           └── mrt_{metric}.sql
├── seeds/                    # Static data
│   └── country_codes.csv
├── snapshots/                # SCD tracking
│   └── snap_{table}.sql
├── target/                   # Compiled SQL (gitignored)
├── tests/
│   ├── generic/
│   └── singular/
├── dbt_project.yml
├── packages.yml
├── profiles.yml              # Connection profiles (gitignored)
├── README.md
└── CLAUDE.md
```

## Airflow Project Structure

### Standard Airflow Layout

```
project/
├── dags/
│   ├── common/               # Shared utilities
│   │   ├── __init__.py
│   │   └── operators.py
│   ├── dag_daily_etl.py
│   ├── dag_weekly_reports.py
│   └── dag_monthly_aggregation.py
├── plugins/                  # Custom plugins
│   └── operators/
├── include/                  # SQL, configs
│   ├── sql/
│   │   └── {domain}/
│   └── configs/
├── tests/
│   └── dags/
│       └── test_dag_daily_etl.py
├── .airflowignore
├── README.md
└── CLAUDE.md
```

## Data Lake Structure

### Medallion Architecture (Bronze/Silver/Gold)

```
data_lake/
├── bronze/                   # Raw ingested data
│   └── {source}/
│       └── {table}/
│           └── year={YYYY}/month={MM}/day={DD}/
│               └── data.parquet
├── silver/                   # Cleaned, conformed
│   └── {domain}/
│       └── {table}/
│           └── year={YYYY}/month={MM}/
│               └── data.parquet
├── gold/                     # Business-ready
│   └── {domain}/
│       └── {table}/
│           └── data.parquet
└── _metadata/
    ├── schemas/
    └── lineage/
```

### Data Lake with Delta Lake

```
data_lake/
├── bronze/
│   └── {table}/
│       ├── _delta_log/       # Transaction log
│       └── part-*.parquet
├── silver/
│   └── {table}/
│       ├── _delta_log/
│       └── part-*.parquet
├── gold/
│   └── {table}/
│       ├── _delta_log/
│       └── part-*.parquet
└── checkpoints/              # Streaming checkpoints
```

## Spark Project Structure

```
project/
├── src/
│   └── {package}/
│       ├── __init__.py
│       ├── jobs/             # Spark jobs
│       │   ├── __init__.py
│       │   └── daily_etl.py
│       ├── transformations/  # Transform functions
│       │   └── __init__.py
│       └── utils/
│           └── __init__.py
├── tests/
│   ├── conftest.py           # Spark session fixture
│   └── jobs/
│       └── test_daily_etl.py
├── configs/
│   └── spark-defaults.conf
├── scripts/
│   └── submit.sh
├── README.md
├── CLAUDE.md
└── pyproject.toml
```

## Directory Purposes

### Data Directories

| Directory | Purpose | Characteristics |
|-----------|---------|-----------------|
| `data/raw/` | Source data | Immutable, partitioned |
| `data/processed/` | Transformed | Clean, conformed |
| `data/features/` | ML features | Versioned, typed |
| `data/external/` | Third-party | Documented source |

### Pipeline Directories

| Directory | Purpose | Contents |
|-----------|---------|----------|
| `pipelines/extract/` | Data extraction | Source connectors |
| `pipelines/transform/` | Transformation | Business logic |
| `pipelines/load/` | Loading | Destination writers |
| `dags/` | Orchestration | Airflow DAGs |

### Configuration

| Directory | Purpose | Contents |
|-----------|---------|----------|
| `configs/` | App config | Environment configs |
| `schemas/` | Data schemas | JSON Schema, Avro |
| `connections/` | Connections | Credentials (gitignored) |

## File Naming Conventions

### Pipeline Files

| Type | Convention | Example |
|------|------------|---------|
| Extract | `extract_{source}.py` | `extract_salesforce.py` |
| Transform | `transform_{domain}.py` | `transform_orders.py` |
| Load | `load_{destination}.py` | `load_snowflake.py` |
| DAG | `dag_{frequency}_{domain}.py` | `dag_daily_sales.py` |

### dbt Models

| Layer | Convention | Example |
|-------|------------|---------|
| Staging | `stg_{source}__{table}.sql` | `stg_shopify__orders.sql` |
| Intermediate | `int_{domain}__{desc}.sql` | `int_orders__aggregated.sql` |
| Marts | `fct_{event}.sql`, `dim_{entity}.sql` | `fct_orders.sql`, `dim_customers.sql` |

### Data Files

| Type | Convention | Example |
|------|------------|---------|
| Raw | `{source}-{YYYYMMDD}-{table}.parquet` | `salesforce-20240115-accounts.parquet` |
| Processed | `{domain}-{table}-{version}.parquet` | `sales-orders-v1.parquet` |
| Partitioned | `{table}/year={YYYY}/month={MM}/` | `orders/year=2024/month=01/` |

## Configuration Management

### Environment Configs

```yaml
# configs/config-prod.yaml
database:
  host: prod-db.example.com
  port: 5432
  name: analytics

warehouse:
  type: snowflake
  account: company.snowflake
  database: ANALYTICS
```

### Secrets Management

Never commit secrets. Use:
- Environment variables
- AWS Secrets Manager
- HashiCorp Vault
- `.env` files (gitignored)

## Testing Strategy

### Test Types

| Type | Location | Purpose |
|------|----------|---------|
| Unit | `tests/unit/` | Transform logic |
| Integration | `tests/integration/` | Pipeline flows |
| Data Quality | `tests/data/` | Schema, values |
| Contract | `tests/contract/` | API contracts |

### Data Quality Tests

```python
# tests/data/test_orders.py
def test_orders_no_nulls():
    """Order ID should never be null"""
    df = spark.read.parquet("data/processed/orders")
    assert df.filter(col("order_id").isNull()).count() == 0

def test_orders_valid_amounts():
    """Order amounts should be positive"""
    df = spark.read.parquet("data/processed/orders")
    assert df.filter(col("amount") <= 0).count() == 0
```

## Anti-Patterns

| Pattern | Problem | Solution |
|---------|---------|----------|
| Hardcoded paths | Environment-specific | Use config files |
| No schema validation | Data drift | Validate schemas |
| Monolithic DAGs | Hard to maintain | Modular DAGs |
| No lineage | Unknown dependencies | Track lineage |
| Secrets in code | Security risk | Use secrets manager |

## Scaling Guidelines

### Small -> Medium

Add:
- Proper configs per environment
- Schema validation
- Data quality tests

### Medium -> Large

Add:
- Orchestration (Airflow)
- Data catalog
- Lineage tracking
- Monitoring/alerting

### Large -> Enterprise

Add:
- Data mesh patterns
- Domain ownership
- Self-service platform
- Governance layer

## References

- [Medallion Architecture](https://databricks.com/glossary/medallion-architecture)
- [dbt Best Practices](https://docs.getdbt.com/guides/best-practices)
- [Airflow Best Practices](https://airflow.apache.org/docs/apache-airflow/stable/best-practices.html)
- [Data Engineering Cookbook](https://github.com/andkret/Cookbook)
