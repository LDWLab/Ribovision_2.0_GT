import loadjs from 'loadjs';
export default (function (deps, id) {
  return new Promise(function (resolve) {
    if (loadjs.isDefined(id)) return resolve();
    loadjs(deps, id, {
      success: resolve,
      error: function error(err) {
        console.error('Deps not found:', err);
        resolve();
      }
    });
  });
});