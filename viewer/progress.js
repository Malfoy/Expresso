// One operation at a time; short operations never show a progress indicator.
export function createProgress(root, bar, label, delay=500) {
  let timer,active=false,done=0,total=0,message='Loading…';
  function render() {
    label.textContent=message;
    if (total>0) {
      bar.max=total;bar.value=Math.min(done,total);
      label.textContent+=` · ${Math.floor(100*Math.min(done,total)/total)}%`;
    } else bar.removeAttribute('value');
  }
  return {
    start() {
      clearTimeout(timer);active=true;done=0;total=0;message='Loading…';root.hidden=true;
      timer=setTimeout(()=>{if (active) {render();root.hidden=false;}},delay);
    },
    plan(count,text) {if (!active) return;total=count;done=0;if (text) message=text;render();},
    advance() {if (active) {done++;render();}},
    message(text) {if (active) {message=text;render();}},
    stop() {active=false;clearTimeout(timer);timer=undefined;root.hidden=true;},
  };
}
