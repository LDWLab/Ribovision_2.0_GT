<template functional>  
  <div class="loading">
    <div class="dots"><span></span> <span></span> <span></span></div>
    <span class="loading-text">Loading Gist</span>
  </div>
</template>

<script>
export default {
  name: 'Loading'
};
</script>

<style scoped>
.loading {
  display: flex;
  align-items: center;
}

.dots {
  pointer-events: none;
  color: transparent;
  font-size: 0
}

.dots span {
    -webkit-animation-name: blink;
            animation-name: blink;
    -webkit-animation-duration: 1.4s;
            animation-duration: 1.4s;
    -webkit-animation-iteration-count: infinite;
            animation-iteration-count: infinite;
    -webkit-animation-fill-mode: both;
            animation-fill-mode: both;
    width: 4px;
    height: 4px;
    background-color: rgb(68, 68, 68);
    display: inline-block;
    border-radius: 50%;
    margin: 0px 1px
  }

.dots span:nth-child(2) {
      -webkit-animation-delay: 0.2s;
              animation-delay: 0.2s;
    }

.dots span:nth-child(3) {
      -webkit-animation-delay: 0.4s;
              animation-delay: 0.4s;
    }

.loading-text {
  margin-left: 8px;
  font-size: 11px;
  opacity: 0.7;
}

@-webkit-keyframes blink {
  0% {
    opacity: 0.2;
  }
  20% {
    opacity: 1;
  }
  100% {
    opacity: 0.2;
  }
}

@keyframes blink {
  0% {
    opacity: 0.2;
  }
  20% {
    opacity: 1;
  }
  100% {
    opacity: 0.2;
  }
}
</style>