#!/usr/bin/env node

'use strict';

// ---------------------------------------------------------------------------
// logpretty — Pretty-print JSON logs in the terminal
// Zero dependencies. Supports pino, winston, and bunyan formats.
// ---------------------------------------------------------------------------

const readline = require('readline');

// ---- ANSI helpers ---------------------------------------------------------

const ANSI = {
  reset:   '\x1b[0m',
  bold:    '\x1b[1m',
  dim:     '\x1b[2m',
  red:     '\x1b[31m',
  green:   '\x1b[32m',
  yellow:  '\x1b[33m',
  blue:    '\x1b[34m',
  magenta: '\x1b[35m',
  cyan:    '\x1b[36m',
  white:   '\x1b[37m',
  gray:    '\x1b[90m',
};

let useColor = true;

function c(code, text) {
  return useColor ? `${code}${text}${ANSI.reset}` : text;
}

// ---- Argument parsing -----------------------------------------------------

function parseArgs(argv) {
  const opts = {
    noColor: false,
    levelKey: null,
    msgKey: null,
    timeKey: null,
    help: false,
  };

  let i = 2; // skip node + script
  while (i < argv.length) {
    const arg = argv[i];
    if (arg === '--no-color') {
      opts.noColor = true;
    } else if (arg === '--level-key' && i + 1 < argv.length) {
      opts.levelKey = argv[++i];
    } else if (arg === '--msg-key' && i + 1 < argv.length) {
      opts.msgKey = argv[++i];
    } else if (arg === '--time-key' && i + 1 < argv.length) {
      opts.timeKey = argv[++i];
    } else if (arg === '--help' || arg === '-h') {
      opts.help = true;
    } else if (arg.startsWith('--level-key=')) {
      opts.levelKey = arg.split('=').slice(1).join('=');
    } else if (arg.startsWith('--msg-key=')) {
      opts.msgKey = arg.split('=').slice(1).join('=');
    } else if (arg.startsWith('--time-key=')) {
      opts.timeKey = arg.split('=').slice(1).join('=');
    }
    i++;
  }
  return opts;
}

function printHelp() {
  const help = `
logpretty — Pretty-print JSON logs in the terminal

USAGE
  cat app.log | logpretty [options]
  docker logs myapp | logpretty
  node server.js | logpretty

OPTIONS
  --no-color          Disable colored output
  --level-key <key>   Custom field name for log level  (default: "level")
  --msg-key <key>     Custom field name for message    (default: "msg", "message")
  --time-key <key>    Custom field name for timestamp  (default: "time", "timestamp", "@timestamp")
  -h, --help          Show this help message

SUPPORTED FORMATS
  pino      {"level":30,"time":1700000000000,"msg":"hello"}
  winston   {"level":"info","message":"hello","timestamp":"..."}
  bunyan    {"level":30,"time":"...","msg":"hello","name":"app"}

EXAMPLES
  cat app.log | logpretty
  cat app.log | logpretty --no-color
  cat app.log | logpretty --level-key severity --msg-key text
`;
  console.log(help.trim());
}

// ---- Level normalization --------------------------------------------------

// Bunyan / pino use numeric levels
const NUMERIC_LEVELS = {
  10: 'TRACE',
  20: 'DEBUG',
  30: 'INFO',
  40: 'WARN',
  50: 'ERROR',
  60: 'FATAL',
};

function normalizeLevel(raw) {
  if (raw == null) return null;
  if (typeof raw === 'number') {
    return NUMERIC_LEVELS[raw] || (raw >= 50 ? 'ERROR' : raw >= 40 ? 'WARN' : raw >= 30 ? 'INFO' : 'DEBUG');
  }
  return String(raw).toUpperCase();
}

function levelColor(level) {
  switch (level) {
    case 'FATAL': return ANSI.red + ANSI.bold;
    case 'ERROR': return ANSI.red;
    case 'WARN':  return ANSI.yellow;
    case 'INFO':  return ANSI.green;
    case 'DEBUG': return ANSI.blue;
    case 'TRACE': return ANSI.gray;
    default:      return ANSI.white;
  }
}

// ---- Timestamp formatting -------------------------------------------------

function formatTimestamp(raw) {
  if (raw == null) return null;
  try {
    // Numeric epoch ms (pino / bunyan)
    if (typeof raw === 'number') {
      return new Date(raw).toISOString().replace('T', ' ').replace('Z', '');
    }
    // ISO string or any Date-parseable string
    const d = new Date(raw);
    if (!isNaN(d.getTime())) {
      return d.toISOString().replace('T', ' ').replace('Z', '');
    }
  } catch (_) { /* fall through */ }
  return String(raw);
}

// ---- Key resolution helpers -----------------------------------------------

function resolveKey(obj, candidates) {
  for (const k of candidates) {
    if (obj[k] !== undefined) return k;
  }
  return null;
}

// ---- Pretty-print a single JSON object ------------------------------------

function formatLine(obj, opts) {
  // Determine keys
  const levelCandidates = opts.levelKey ? [opts.levelKey] : ['level'];
  const msgCandidates   = opts.msgKey   ? [opts.msgKey]   : ['msg', 'message'];
  const timeCandidates  = opts.timeKey  ? [opts.timeKey]  : ['time', 'timestamp', '@timestamp'];

  const levelKey = resolveKey(obj, levelCandidates);
  const msgKey   = resolveKey(obj, msgCandidates);
  const timeKey  = resolveKey(obj, timeCandidates);

  const rawLevel = levelKey ? obj[levelKey] : null;
  const level    = normalizeLevel(rawLevel);
  const message  = msgKey ? String(obj[msgKey]) : '';
  const time     = timeKey ? formatTimestamp(obj[timeKey]) : null;

  // Build the main line
  const parts = [];

  if (time) {
    parts.push(c(ANSI.gray, time));
  }

  if (level) {
    const padded = level.padEnd(5);
    parts.push(c(levelColor(level), padded));
  }

  if (message) {
    parts.push(c(ANSI.white + ANSI.bold, message));
  }

  const mainLine = parts.join(' ');

  // Collect extra fields
  const skipKeys = new Set();
  if (levelKey) skipKeys.add(levelKey);
  if (msgKey) skipKeys.add(msgKey);
  if (timeKey) skipKeys.add(timeKey);
  // Also skip common meta keys that are noise
  const metaKeys = ['pid', 'hostname', 'name', 'v'];
  for (const mk of metaKeys) {
    if (obj[mk] !== undefined) skipKeys.add(mk);
  }

  const extraKeys = Object.keys(obj).filter((k) => !skipKeys.has(k));

  const lines = [mainLine];

  if (extraKeys.length > 0) {
    for (const k of extraKeys) {
      let val = obj[k];
      if (typeof val === 'object' && val !== null) {
        val = JSON.stringify(val, null, 2)
          .split('\n')
          .map((l, i) => (i === 0 ? l : '              ' + l))
          .join('\n');
      }
      lines.push(`    ${c(ANSI.cyan, k)}${c(ANSI.gray, ':')} ${val}`);
    }
  }

  return lines.join('\n');
}

// ---- Main -----------------------------------------------------------------

function main() {
  const opts = parseArgs(process.argv);

  if (opts.help) {
    printHelp();
    process.exit(0);
  }

  if (opts.noColor || process.env.NO_COLOR !== undefined) {
    useColor = false;
  }

  // Check if stdin is a TTY (no pipe)
  if (process.stdin.isTTY) {
    console.error('logpretty: no input detected. Pipe JSON logs into logpretty.');
    console.error('');
    console.error('  Example:  cat app.log | logpretty');
    console.error('  Help:     logpretty --help');
    process.exit(1);
  }

  const rl = readline.createInterface({
    input: process.stdin,
    crlfDelay: Infinity,
  });

  rl.on('line', (line) => {
    const trimmed = line.trim();
    if (!trimmed) {
      console.log(line);
      return;
    }

    // Attempt JSON parse
    let obj;
    try {
      obj = JSON.parse(trimmed);
    } catch (_) {
      // Not JSON — pass through unchanged
      console.log(line);
      return;
    }

    if (typeof obj !== 'object' || obj === null || Array.isArray(obj)) {
      // Valid JSON but not an object — pass through
      console.log(line);
      return;
    }

    try {
      console.log(formatLine(obj, opts));
    } catch (_) {
      console.log(line);
    }
  });

  rl.on('close', () => {
    process.exit(0);
  });

  // Handle broken pipe gracefully
  process.stdout.on('error', (err) => {
    if (err.code === 'EPIPE') {
      process.exit(0);
    }
    throw err;
  });
}

main();
