export var isExternalLink = function isExternalLink(link) {
  return /^https?:\/\//.test(link);
};
export var slugify = function slugify(str) {
  var RE = /[\u2000-\u206F\u2E00-\u2E7F\\'!"#$%&()*+,./:;<=>?@[\]^`{|}~]/g;
  var REPLACEMENT = '-';
  var WHITESPACE = /\s/g;
  return str.trim().replace(RE, '').replace(WHITESPACE, REPLACEMENT).toLowerCase();
};
export var getFileUrl = function getFileUrl(sourcePath, path) {
  sourcePath = sourcePath || '.'; // Remove trailing slash in `sourcePath`
  // Since `path` always starts with slash

  sourcePath = sourcePath.replace(/\/$/, '');
  var result = sourcePath + path;
  return result.replace(/^\.\//, '');
};
export var getFilenameByPath = function getFilenameByPath(path) {
  // Ensure path always starts with slash
  path = path.replace(/^\/?/, '/'); // Add .md suffix

  if (!/\.md$/.test(path)) {
    path = /\/$/.test(path) ? path + "README.md" : path + ".md";
  }

  return path;
};
export var inBrowser = typeof window !== 'undefined';