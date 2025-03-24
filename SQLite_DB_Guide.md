## Direct SQLite Database Access in Docker Container

Execute these commands on your droplet to directly access the SQLite database in the running container:

```bash
# Enter the running backend container with bash
docker exec -it optimol-backend bash

# Locate the database file - it's likely in /app/instance
find /app -name "*.db"

# Install SQLite if not already in container (only needed once)
apt-get update && apt-get install -y sqlite3

# Access the SQLite database (adjust path if find returns different location)
sqlite3 /app/instance/users.db
```

## Technical SQLite Commands

Once in the SQLite prompt:

```sql
-- Show SQLite version
.version

-- Configure output for readability
.mode column
.headers on
.width 36 10 15 0

-- List all tables
.tables

-- Inspect schemas
.schema optimization
.schema user

-- Get column definitions
PRAGMA table_info(optimization);
PRAGMA table_info(user);

-- View records
SELECT id, user_id, optimization_type, created_at FROM optimization ORDER BY created_at DESC LIMIT 10;
SELECT id, email, subscription_status FROM user LIMIT 10;

-- Count records in tables
SELECT COUNT(*) FROM optimization;
SELECT COUNT(*) FROM user;

-- View foreign keys
PRAGMA foreign_key_list(optimization);

-- Execute schema modifications
PRAGMA foreign_keys=off;
BEGIN TRANSACTION;
CREATE TABLE optimization_new (
  id VARCHAR(36) PRIMARY KEY,
  user_id INTEGER REFERENCES user(id),
  optimization_type VARCHAR(50) NOT NULL,
  parameters TEXT NOT NULL,
  result TEXT,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
INSERT INTO optimization_new SELECT * FROM optimization;
DROP TABLE optimization;
ALTER TABLE optimization_new RENAME TO optimization;
COMMIT;
PRAGMA foreign_keys=on;

-- Data manipulation operations
DELETE FROM optimization WHERE id = 'specific-uuid';
UPDATE user SET subscription_status = 'active' WHERE id = 1;

## Database Backup and Restore

# Inside the temporary container, create a backup
sqlite3 /data/users.db .dump > /data/backup.sql

# Restore from backup (caution: overwrites existing database)
cat /data/backup.sql | sqlite3 /data/users.db.new


-- Exit SQLite
.quit
```

Exit the container shell when finished:
```bash
exit
```





