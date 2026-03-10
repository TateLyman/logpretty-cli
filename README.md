# logpretty

![Version](https://img.shields.io/badge/version-1.0.0-blue) ![License](https://img.shields.io/badge/license-MIT-green)

> Pretty-print JSON logs in the terminal. Supports pino, winston, and bunyan formats out of the box.

If you use structured logging, your raw logs look like this:

```
{"level":30,"time":1700000000000,"pid":12345,"hostname":"server1","msg":"User signed in","userId":"u_8fk2","latency":42}
{"level":50,"time":1700000001000,"pid":12345,"hostname":"server1","msg":"Database connection failed","error":"ECONNREFUSED","retries":3}
{"level":40,"time":1700000002000,"pid":12345,"hostname":"server1","msg":"Rate limit approaching","current":980,"max":1000}
```

Pipe them through `logpretty` and get this:

```
2023-11-14 22:13:20.000 INFO  User signed in
    userId: u_8fk2
    latency: 42
2023-11-14 22:13:21.000 ERROR Database connection failed
    error: ECONNREFUSED
    retries: 3
2023-11-14 22:13:22.000 WARN  Rate limit approaching
    current: 980
    max: 1000
```

Timestamps in gray, levels color-coded (ERROR = red, WARN = yellow, INFO = green, DEBUG = blue), messages in bold white, and extra fields neatly indented below.

## Installation

```bash
# npm
npm install -g logpretty-cli

# yarn
yarn global add logpretty-cli

# pnpm
pnpm add -g logpretty-cli
```

## Usage

Pipe any JSON log stream through `logpretty`:

```bash
# From a log file
cat app.log | logpretty

# From a running process
node server.js 2>&1 | logpretty

# Docker logs
docker logs -f myapp | logpretty

# PM2 logs
pm2 logs --raw | logpretty
```

Non-JSON lines pass through unchanged, so mixed output works fine.

## Options

```
--no-color          Disable colored output
--level-key <key>   Custom field name for log level  (default: "level")
--msg-key <key>     Custom field name for message    (default: "msg", "message")
--time-key <key>    Custom field name for timestamp  (default: "time", "timestamp", "@timestamp")
-h, --help          Show help
```

## Supported Log Formats

### pino

```json
{"level":30,"time":1700000000000,"pid":42,"hostname":"app","msg":"Request completed","responseTime":12}
```

### winston (JSON transport)

```json
{"level":"info","message":"Server started","timestamp":"2023-11-14T22:13:20.000Z","port":3000}
```

### bunyan

```json
{"name":"myapp","hostname":"server","pid":1234,"level":30,"msg":"Listening","time":"2023-11-14T22:13:20.000Z","v":0,"port":8080}
```

### Custom formats

Use the `--level-key`, `--msg-key`, and `--time-key` flags for any format:

```bash
cat custom.log | logpretty --level-key severity --msg-key text --time-key ts
```

## How It Works

1. Reads stdin line by line
2. Attempts to parse each line as JSON
3. If valid JSON object: extracts timestamp, level, and message fields, then pretty-prints with ANSI colors. Extra fields are displayed indented below the main line.
4. If not JSON: prints the line unchanged

## Zero Dependencies

This tool has **no dependencies**. It uses only Node.js built-in modules (`readline`). Install it and it just works.

## License

[MIT](./LICENSE) -- Tate Lyman
